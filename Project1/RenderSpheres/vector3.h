#ifndef VECTOR3_H
#define VECTOR3_H

class Vector3
{
private:
public:
	double v[3];
	Vector3(double v0, double v1, double v2) : v{v0, v1, v2} {};
	Vector3() : v{0, 0, 0} {};

	Vector3 operator-() const { return Vector3(-v[0], -v[1], -v[2]); }
	Vector3 &operator+=(const Vector3 u)
	{
		v[0] += u.v[0];
		v[1] += u.v[1];
		v[2] += u.v[2];
		return *this;
	}

	Vector3 &operator*=(double k)
	{
		v[0] *= k;
		v[1] *= k;
		v[2] *= k;
		return *this;
	}

	Vector3 &operator/=(double k)
	{
		return *this *= 1 / k;
	}
};

Vector3 operator+(const Vector3 &u, const Vector3 &w)
{
	return Vector3(u.v[0] + w.v[0], u.v[1] + w.v[1], u.v[2] + w.v[2]);
}

double dot(const Vector3 &u, const Vector3 &w)
{
	return u.v[0] * w.v[0] + u.v[1] * w.v[1] + u.v[2] * w.v[2];
}

Vector3 cross(const Vector3 &u, const Vector3 &w)
{
	return Vector3(u.v[1] * w.v[2] - u.v[2] * w.v[1],
				   u.v[2] * w.v[0] - u.v[0] * w.v[2],
				   u.v[0] * w.v[1] - u.v[1] * w.v[0]);
}

#endif