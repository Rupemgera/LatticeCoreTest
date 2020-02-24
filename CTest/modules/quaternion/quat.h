#pragma once
#include <eigen3/Eigen/Dense>

typedef Eigen::Matrix<double, 4, 1> Vector4f;
typedef Eigen::Matrix<double, 3, 1> Vector3f;
typedef Eigen::Matrix<double, 3, 3> Matrix3f;

class Quaternion : public Eigen::Matrix<double, 4, 1> {
public:
  typedef Eigen::Matrix<double, 4, 1> Base;

  Quaternion() : Base() {}
  Quaternion(double x, double y, double z, double w) : Base(x, y, z, w) {}
  Quaternion(const Eigen::Matrix<double, 3, 1> &xyz, double w)
      : Base(xyz.x(), xyz.y(), xyz.z(), w) {}

  /// Construct a vector from MatrixBase (needed to play nice with Eigen)
  template <typename Derived>
  Quaternion(const Eigen::MatrixBase<Derived> &p) : Base(p) {}

  /// Assign a vector from MatrixBase (needed to play nice with Eigen)
  template <typename Derived>
  Quaternion &operator=(const Eigen::MatrixBase<Derived> &p) {
    this->Base::operator=(p);
    return *this;
  }

  friend Quaternion operator*(const Quaternion &q0, const Quaternion &q1) {
    Quaternion result;
    result << q0.head<3>().cross(q1.head<3>()) + q1.w() * q0.head<3>() +
                  q0.w() * q1.head<3>(),
        q0.w() * q1.w() - q0.head<3>().dot(q1.head<3>());
    return result;
  }

  Quaternion complement() const { return Quaternion(-x(), -y(), -z(), w()); }

  Quaternion exp() const {
    double ri = head<3>().norm();
    double sinRi = std::sin(ri);
    double cosRi = std::cos(ri);
    double expW = std::exp(w());

    return (Quaternion() << head<3>() * (expW * sinRi / ri), expW * cosRi)
        .finished();
  }

  Quaternion log() const {
    double ri = head<3>().norm();
    double rq = norm();
    double acosR = std::acos(w() / rq);
    double logRq = std::log(rq);
    return (Quaternion() << head<3>() * (acosR / ri), logRq).finished();
  }

  Vector4f pow(double exponent) const {
    return Quaternion(log() * exponent).exp();
  }

  Eigen::Matrix<double, 3, 3> toMatrix() const {
    double xx = x() * x(), yy = y() * y(), zz = z() * z();
    double xy = x() * y(), xz = x() * z(), yz = y() * z();
    double xw = x() * w(), yw = y() * w(), zw = z() * w();

    return (Matrix3f() << 1 - 2 * yy - 2 * zz, 2 * (xy - zw), 2 * (xz + yw),
            2 * (xy + zw), 1 - 2 * xx - 2 * zz, 2 * (yz - xw), 2 * (xz - yw),
            2 * (yz + xw), 1 - 2 * xx - 2 * yy)
        .finished();
  }

  Quaternion align(const Vector3f &n_) const {
    Quaternion n = complement() * (Quaternion(n_, 0.f) * *this);
    Quaternion na = n.cwiseAbs();

    Quaternion s;
    if (na.x() > na.y() && na.x() > na.z())
      s = Quaternion(0, -n.z(), n.y(), n.x());
    else if (na.y() > na.z())
      s = Quaternion(n.z(), 0.f, -n.x(), n.y());
    else
      s = Quaternion(-n.y(), n.x(), 0.f, n.z());
    s.w() += std::copysign(1.f, s.w());
    s.normalize();

    return *this * s;
  }

  static Quaternion findRotation(const Quaternion &q, const Quaternion &ref) {
    double dp[4] = {ref.dot(Quaternion(q.w(), q.z(), -q.y(), -q.x())),
                    ref.dot(Quaternion(-q.z(), q.w(), q.x(), -q.y())),
                    ref.dot(Quaternion(q.y(), -q.x(), q.w(), -q.z())),
                    ref.dot(Quaternion(q.x(), q.y(), q.z(), q.w()))};

    double a[4] = {std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]),
                   std::abs(dp[3])};

    // find M, N, max indices of a[i]
    int M = 0, N = 1, m = 2, n = 3;
    if (a[M] < a[m])
      std::swap(M, m);
    if (a[N] < a[n])
      std::swap(N, n);
    if (a[M] < a[N]) {
      std::swap(M, N);
      m = n;
    }
    if (a[N] < a[m])
      N = m;

    const double s = std::sqrt(.5f);
    const double h = .5f;

    double vA = a[M];
    double vB = (a[M] + a[N]) * s;
    double vC = (a[0] + a[1] + a[2] + a[3]) * h;

    Quaternion result = Quaternion::Zero();
    if (vA > vB && vA > vC) {
      result[M] = std::copysign(1.f, dp[M]);
    } else if (vB > vC) {
      result[M] = std::copysign(s, dp[M]);
      result[N] = std::copysign(s, dp[N]);
    } else {
      result[0] = std::copysign(h, dp[0]);
      result[1] = std::copysign(h, dp[1]);
      result[2] = std::copysign(h, dp[2]);
      result[3] = std::copysign(h, dp[3]);
    }

    return result;
  }

  static Quaternion applyRotation(const Quaternion &q, const Quaternion &ref) {
    Quaternion qi[4] = {Quaternion(q.w(), q.z(), -q.y(), -q.x()),
                        Quaternion(-q.z(), q.w(), q.x(), -q.y()),
                        Quaternion(q.y(), -q.x(), q.w(), -q.z()),
                        Quaternion(q.x(), q.y(), q.z(), q.w())};
    double dp[4] = {ref.dot(qi[0]), ref.dot(qi[1]), ref.dot(qi[2]),
                    ref.dot(qi[3])};
    double a[4] = {std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]),
                   std::abs(dp[3])};

    // find M, N, max indices of a[i]
    int M = 0, N = 1, m = 2, n = 3;
    if (a[M] < a[m])
      std::swap(M, m);
    if (a[N] < a[n])
      std::swap(N, n);
    if (a[M] < a[N]) {
      std::swap(M, N);
      m = n;
    }
    if (a[N] < a[m])
      N = m;

    const double s = std::sqrt(.5f);
    const double h = .5f;

    double vA = a[M];
    double vB = (a[M] + a[N]) * s;
    double vC = (a[0] + a[1] + a[2] + a[3]) * h;

    if (vA > vB && vA > vC) {
      return qi[M] * std::copysign(1.f, dp[M]);
    } else if (vB > vC) {
      return qi[M] * std::copysign(s, dp[M]) + qi[N] * std::copysign(s, dp[N]);
    } else {
      return qi[0] * std::copysign(h, dp[0]) + qi[1] * std::copysign(h, dp[1]) +
             qi[2] * std::copysign(h, dp[2]) + qi[3] * std::copysign(h, dp[3]);
    }
  }
};
