#include <Novice.h>
#include <cmath>
#include <cassert>
#include <imgui.h>
#define _USE_MATH_DEFINES
#include "math.h"

const char kWindowTitle[] = "LE2B_12_サクライショウセイ_タイトル";

struct Vector3 {
	float x, y, z;
};

struct Matrix4x4 {
	float m[4][4];
};

struct Plane {
	Vector3 normal; //!< 法線
	float distance; //!< 距離
};

struct Segment {
	Vector3 origin;
	Vector3 diff;
};

struct Triangle{
	Vector3 vertices[3];
};

Vector3 Add(const Vector3& v1, const Vector3& v2);
Vector3 Subtract(const Vector3& v1, const Vector3& v2);
float Dot(const Vector3& v1, const Vector3& v2);
Vector3 Normalize(const Vector3& v);
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2);
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate);
Matrix4x4 Inverse(const Matrix4x4& m);
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip);
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth);
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix);
Vector3 Cross(const Vector3& v1, const Vector3& v2);
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix);
void DrawSegment(const Segment& segment, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color);
bool IsCollision(const Segment& segment, const Plane& plane);
void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color);

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };

	Vector3 cameraRotate{ 0.26F,0.0F,0.0F };

	uint32_t segmentcolor = WHITE;;

	Segment segment = {};

	Triangle triangle{};

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("Segment.origin", &segment.origin.x, 0.01f);
		ImGui::DragFloat3("Segment.diff", &segment.diff.x, 0.01f);
		ImGui::DragFloat3("vertices[0]", &triangle.vertices[0].x, 0.01f);
		ImGui::DragFloat3("vertices[1]", &triangle.vertices[1].x, 0.01f);
		ImGui::DragFloat3("vertices[2]", &triangle.vertices[2].x, 0.01f);
		ImGui::End();

		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f, 1.0f,1.0f }, cameraRotate, cameraTranslate);

		Matrix4x4 viewMatrix = Inverse(cameraMatrix);

		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, 1280.0f / 720.0f, 0.1f, 100.f);

		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, 1280.0f, 720.0f, 0.0f, 1.0f);

		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);

		Vector3 v1 = Subtract(triangle.vertices[1], triangle.vertices[0]);
		Vector3 v2 = Subtract(triangle.vertices[2], triangle.vertices[1]);

		Vector3 v = Cross(v1, v2);

		v = Normalize(v);

		float d = Dot(triangle.vertices[0], v);

		Plane trianglePlane = { v,d };

		float dot = Dot(trianglePlane.normal, segment.diff);

		float t = (trianglePlane.distance - Dot(segment.origin, trianglePlane.normal)) / dot;

		Vector3 p = Add(segment.origin, { segment.diff.x * t,segment.diff.y * t,segment.diff.z * t });

		if (IsCollision(segment, trianglePlane))
		{
			Vector3 v01 = Subtract(triangle.vertices[0], triangle.vertices[1]);
			Vector3 v12 = Subtract(triangle.vertices[1], triangle.vertices[2]);
			Vector3 v20 = Subtract(triangle.vertices[2], triangle.vertices[0]);

			Vector3 vp0 = Subtract(triangle.vertices[0], p);
			Vector3 vp1 = Subtract(triangle.vertices[1], p);
			Vector3 vp2 = Subtract(triangle.vertices[2], p);

			Vector3 cross01 = Cross(v01, vp1);
			Vector3 cross12 = Cross(v12, vp2);
			Vector3 cross20 = Cross(v20, vp0);

			if (Dot(cross01, trianglePlane.normal) >= 0.0f &&
				Dot(cross12, trianglePlane.normal) >= 0.0f &&
				Dot(cross20, trianglePlane.normal) >= 0.0f) {
				segmentcolor = RED;
			} else {
				segmentcolor = WHITE;
			}
		} else {
			segmentcolor = WHITE;
		}

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(viewProjectionMatrix, viewportMatrix);

		DrawSegment(segment, viewProjectionMatrix, viewportMatrix, segmentcolor);

		DrawTriangle(triangle, viewProjectionMatrix, viewportMatrix, WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}

Vector3 Add(const Vector3& v1, const Vector3& v2)
{
	Vector3 a;
	a.x = v1.x + v2.x;
	a.y = v1.y + v2.y;
	a.z = v1.z + v2.z;
	return a;
}

Vector3 Subtract(const Vector3& v1, const Vector3& v2)
{
	Vector3 a;
	a.x = v1.x - v2.x;
	a.y = v1.y - v2.y;
	a.z = v1.z - v2.z;
	return a;
}

float Dot(const Vector3& v1, const Vector3& v2)
{
	float b;
	b = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
	return b;
}

Vector3 Normalize(const Vector3& v)
{
	Vector3 a;
	float b;
	b = sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
	a.x = v.x / b;
	a.y = v.y / b;
	a.z = v.z / b;
	return a;
}

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 a;
	a.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
	a.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
	a.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
	a.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

	a.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
	a.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
	a.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
	a.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

	a.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
	a.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
	a.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
	a.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

	a.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
	a.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
	a.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
	a.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];
	return a;
}

Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 a;
	a.m[0][0] = 1;
	a.m[0][1] = 0;
	a.m[0][2] = 0;
	a.m[0][3] = 0;
	a.m[1][0] = 0;
	a.m[1][1] = std::cos(radian);
	a.m[1][2] = std::sin(radian);
	a.m[1][3] = 0;
	a.m[2][0] = 0;
	a.m[2][1] = -std::sin(radian);
	a.m[2][2] = std::cos(radian);
	a.m[2][3] = 0;
	a.m[3][0] = 0;
	a.m[3][1] = 0;
	a.m[3][2] = 0;
	a.m[3][3] = 1;
	return a;
}

Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 a;
	a.m[0][0] = std::cos(radian);
	a.m[0][1] = 0;
	a.m[0][2] = -std::sin(radian);
	a.m[0][3] = 0;
	a.m[1][0] = 0;
	a.m[1][1] = 1;
	a.m[1][2] = 0;
	a.m[1][3] = 0;
	a.m[2][0] = std::sin(radian);
	a.m[2][1] = 0;
	a.m[2][2] = std::cos(radian);
	a.m[2][3] = 0;
	a.m[3][0] = 0;
	a.m[3][1] = 0;
	a.m[3][2] = 0;
	a.m[3][3] = 1;
	return a;
}

Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 a;
	a.m[0][0] = std::cos(radian);
	a.m[0][1] = std::sin(radian);
	a.m[0][2] = 0;
	a.m[0][3] = 0;
	a.m[1][0] = -std::sin(radian);
	a.m[1][1] = std::cos(radian);
	a.m[1][2] = 0;
	a.m[1][3] = 0;
	a.m[2][0] = 0;
	a.m[2][1] = 0;
	a.m[2][2] = 1;
	a.m[2][3] = 0;
	a.m[3][0] = 0;
	a.m[3][1] = 0;
	a.m[3][2] = 0;
	a.m[3][3] = 1;
	return a;
}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 a;
	Matrix4x4 b;

	b = Multiply(MakeRotateXMatrix(rotate.x), Multiply(MakeRotateYMatrix(rotate.y), MakeRotateZMatrix(rotate.z)));
	a.m[0][0] = scale.x * b.m[0][0];
	a.m[0][1] = scale.x * b.m[0][1];
	a.m[0][2] = scale.x * b.m[0][2];
	a.m[0][3] = 0;
	a.m[1][0] = scale.y * b.m[1][0];
	a.m[1][1] = scale.y * b.m[1][1];
	a.m[1][2] = scale.y * b.m[1][2];
	a.m[1][3] = 0;
	a.m[2][0] = scale.z * b.m[2][0];
	a.m[2][1] = scale.z * b.m[2][1];
	a.m[2][2] = scale.z * b.m[2][2];
	a.m[2][3] = 0;
	a.m[3][0] = translate.x;
	a.m[3][1] = translate.y;
	a.m[3][2] = translate.z;
	a.m[3][3] = 1;

	return a;
}

Matrix4x4 Inverse(const Matrix4x4& m)
{
	Matrix4x4 a;
	float detA =
		m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2] -

		m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2] -

		m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2] +

		m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2] +

		m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2] -

		m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2] -

		m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0] +

		m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	a.m[0][0] = 1 / detA * (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]);
	a.m[0][1] = 1 / detA * (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]);
	a.m[0][2] = 1 / detA * (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]);
	a.m[0][3] = 1 / detA * (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]);

	a.m[1][0] = 1 / detA * -(m.m[1][0] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][2]);
	a.m[1][1] = 1 / detA * -(-m.m[0][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][2]);
	a.m[1][2] = 1 / detA * -(m.m[0][0] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][2]);
	a.m[1][3] = 1 / detA * -(-m.m[0][0] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][2]);

	a.m[2][0] = 1 / detA * (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]);
	a.m[2][1] = 1 / detA * (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]);
	a.m[2][2] = 1 / detA * (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]);
	a.m[2][3] = 1 / detA * (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]);

	a.m[3][0] = 1 / detA * -(m.m[1][0] * m.m[2][1] * m.m[3][2] + m.m[1][1] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][2] - m.m[1][0] * m.m[2][2] * m.m[3][1]);
	a.m[3][1] = 1 / detA * -(-m.m[0][0] * m.m[2][1] * m.m[3][2] - m.m[0][1] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][2] + m.m[0][0] * m.m[2][2] * m.m[3][1]);
	a.m[3][2] = 1 / detA * -(m.m[0][0] * m.m[1][1] * m.m[3][2] + m.m[0][1] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][2] - m.m[0][0] * m.m[1][2] * m.m[3][1]);
	a.m[3][3] = 1 / detA * -(-m.m[0][0] * m.m[1][1] * m.m[2][2] - m.m[0][1] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][2] + m.m[0][0] * m.m[1][2] * m.m[2][1]);

	return a;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
	Matrix4x4 a;
	a.m[0][0] = 1 / aspectRatio * 1 / std::tan(fovY / 2);
	a.m[0][1] = 0;
	a.m[0][2] = 0;
	a.m[0][3] = 0;
	a.m[1][0] = 0;
	a.m[1][1] = 1 / std::tan(fovY / 2);
	a.m[1][2] = 0;
	a.m[1][3] = 0;
	a.m[2][0] = 0;
	a.m[2][1] = 0;
	a.m[2][2] = farClip / (farClip - nearClip);
	a.m[2][3] = 1;
	a.m[3][0] = 0;
	a.m[3][1] = 0;
	a.m[3][2] = -nearClip * farClip / (farClip - nearClip);
	a.m[3][3] = 0;
	return a;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
	Matrix4x4 a;
	a.m[0][0] = width / 2;
	a.m[0][1] = 0;
	a.m[0][2] = 0;
	a.m[0][3] = 0;
	a.m[1][0] = 0;
	a.m[1][1] = -height / 2;
	a.m[1][2] = 0;
	a.m[1][3] = 0;
	a.m[2][0] = 0;
	a.m[2][1] = 0;
	a.m[2][2] = maxDepth - minDepth;
	a.m[2][3] = 0;
	a.m[3][0] = left + width / 2;
	a.m[3][1] = top + height / 2;
	a.m[3][2] = minDepth;
	a.m[3][3] = 1;
	return a;
}

Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix)
{
	Vector3 b;
	b.x = (vector.x * matrix.m[0][0]) + (vector.y * matrix.m[1][0]) + (vector.z * matrix.m[2][0]) + (1.0f * matrix.m[3][0]);
	b.y = (vector.x * matrix.m[0][1]) + (vector.y * matrix.m[1][1]) + (vector.z * matrix.m[2][1]) + (1.0f * matrix.m[3][1]);
	b.z = (vector.x * matrix.m[0][2]) + (vector.y * matrix.m[1][2]) + (vector.z * matrix.m[2][2]) + (1.0f * matrix.m[3][2]);
	float w = (vector.x * matrix.m[0][3]) + (vector.y * matrix.m[1][3]) + (vector.z * matrix.m[2][3]) + (1.0f * matrix.m[3][3]);
	assert(w != 0.0f);
	b.x /= w;
	b.y /= w;
	b.z /= w;
	return b;
}

Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
	Vector3 a;
	a.x = v1.y * v2.z - v1.z * v2.y;
	a.y = v1.z * v2.x - v1.x * v2.z;
	a.z = v1.x * v2.y - v1.y * v2.x;
	return a;
}

bool IsCollision(const Segment& segment, const Plane& plane)
{
	float a = Dot(plane.normal, segment.diff);

	if (a == 0.0f)
	{
		return false;
	}

	float t = (plane.distance - Dot(segment.origin, plane.normal)) / a;

	if (0 <= t && t <= 1)
	{
		return true;
	} else {
		return false;
	}
}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);
	// 奥から手前への線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {

		Matrix4x4 startWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, { (xIndex * kGridEvery),0.0f,kGridHalfWidth });
		Matrix4x4 endWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, { (xIndex * kGridEvery),0.0f,-kGridHalfWidth });

		Matrix4x4 startWorldViewProjectionMatrix = Multiply(startWorldMatrix, viewProjectionMatrix);
		Matrix4x4 endWorldViewProjectionMatrix = Multiply(endWorldMatrix, viewProjectionMatrix);

		Vector3 ndcStartVertex = Transform({ -2,0,0 }, startWorldViewProjectionMatrix);
		Vector3 ndcEndVertex = Transform({ -2,0,0 }, endWorldViewProjectionMatrix);

		Vector3 screenStartVertex = Transform(ndcStartVertex, viewportMatrix);
		Vector3 screenEndVertex = Transform(ndcEndVertex, viewportMatrix);
		Novice::DrawLine(int(screenStartVertex.x), int(screenStartVertex.y), int(screenEndVertex.x), int(screenEndVertex.y), 0xAAAAAAFF);
	}

	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {

		Matrix4x4 startWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, { kGridHalfWidth,0.0f,(zIndex * kGridEvery) });
		Matrix4x4 endWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, { -kGridHalfWidth,0.0f,(zIndex * kGridEvery) });

		Matrix4x4 startWorldViewProjectionMatrix = Multiply(startWorldMatrix, viewProjectionMatrix);
		Matrix4x4 endWorldViewProjectionMatrix = Multiply(endWorldMatrix, viewProjectionMatrix);

		Vector3 ndcStartVertex = Transform({ 0,0,-2 }, startWorldViewProjectionMatrix);
		Vector3 ndcEndVertex = Transform({ 0,0,-2 }, endWorldViewProjectionMatrix);

		Vector3 screenStartVertex = Transform(ndcStartVertex, viewportMatrix);
		Vector3 screenEndVertex = Transform(ndcEndVertex, viewportMatrix);

		Novice::DrawLine(int(screenStartVertex.x), int(screenStartVertex.y), int(screenEndVertex.x), int(screenEndVertex.y), 0xAAAAAAFF);
	}
}

void DrawSegment(const Segment& segment, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 origin = Transform(Transform(segment.origin, viewProjectionMatrix), viewportMatrix);
	Vector3 end = Transform(Transform(Add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);

	Novice::DrawLine(static_cast<int>(origin.x), static_cast<int>(origin.y), static_cast<int>(end.x), static_cast<int>(end.y), color);
}

void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	Vector3 Vertices[3];
	for (int i = 0; i < 3; i++)
	{
		Vertices[i] = Transform(Transform(triangle.vertices[i], viewProjectionMatrix), viewportMatrix);
	}

	Novice::DrawTriangle(static_cast<int>(Vertices[0].x), static_cast<int>(Vertices[0].y),
		static_cast<int>(Vertices[1].x), static_cast<int>(Vertices[1].y),
		static_cast<int>(Vertices[2].x), static_cast<int>(Vertices[2].y), color, kFillModeWireFrame);
}
