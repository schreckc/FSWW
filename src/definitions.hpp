#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <Eigen/Core>

using namespace Eigen;

//#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION

#define FLOAT double
#define VEC2 Vector2d
#define VEC3 Vector3d
#define VEC4 Vector4d
#define VECX VectorXd
#define MAT2 Matrix2d
#define MAT3 Matrix3d
#define MAT4 Matrix4d
#define MATX MatrixXd
#define ANGLE_AXIS AngleAxisd
#define QUATERNION Quaterniond
#define COMPLEX std::complex<double>
#define VEC2C Vector2cd
#define VECXC VectorXcd
#define MATXC MatrixXcd

#else

#define FLOAT float
#define VEC2 Vector2f
#define VEC3 Vector3f
#define VEC4 Vector4f
#define VECX VectorXf
#define MAT2 Matrix2f
#define MAT3 Matrix3f
#define MAT4 Matrix4f
#define MATX MatrixXf
#define ANGLE_AXIS AngleAxisf
#define QUATERNION Quaternionf
#define COMPLEX std::complex<float>
#define VEC2C Vector2cf
#define VECXC VectorXcf
#define MATXC MatrixXcf

#endif

namespace definitions {
  extern COMPLEX i_;
  extern COMPLEX r_;
};


#endif
