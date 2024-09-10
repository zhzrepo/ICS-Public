#ifndef VECTOR3D_HH
#define VECTOR3D_HH

#include <cmath>

class Vector3D {
public:
    double x, y, z;

    // Constructor
    Vector3D(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

    // Vector addition
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    // Vector subtraction
    Vector3D operator-(const Vector3D& v) const {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    // Scalar multiplication: Vector * Scalar
    Vector3D operator*(G4double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    // Scalar multiplication: Scalar * Vector
    friend Vector3D operator*(G4double scalar, const Vector3D& vec) {
        return Vector3D(vec.x * scalar, vec.y * scalar, vec.z * scalar);
    }

    // Scalar division
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }

    // Dot product
    double dot(const Vector3D& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // Cross product
    Vector3D cross(const Vector3D& v) const {
        return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Normalize the vector
    Vector3D normalize() const {
        G4double len = sqrt(x * x + y * y + z * z);
        if (len == 0) return Vector3D(0, 0, 0);  // To handle division by zero if vector is zero
        return Vector3D(x / len, y / len, z / len);
    }

    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const Vector3D& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};

#endif // VECTOR3D_HH
