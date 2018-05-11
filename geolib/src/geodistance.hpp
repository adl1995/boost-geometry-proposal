/**
 * @file geodistance.hpp
 * @author Adeel Ahmad
 *
 * Geographical distance between two points on the Earth's surface.
 * The current implementation includes the Haversine formula.
 *
 */
#ifndef GEOLIB_SRC_GEODISTANCE_HPP
#define GEOLIB_SRC_GEODISTANCE_HPP

#include "geographiclib-geodesic.hpp"

// Define a traits system which extends a generic getRadian
// and getDegree function to be used for accessing its
// elements using the AccessPoint structure.
namespace PointTrait
{
  template <typename P, int D>
  struct AccessPoint {};

  // These free functions enables us to call getRadian<0>(point),
  // where 0 represents latitude. This works for any point
  // having the AccessPoint specialization.
  template <int D, typename P>
  inline double getRadian(P const& p)
  {
      return AccessPoint<P, D>::getRadian(p);
  }

  template <int D, typename P>
  inline double getDegree(P const& p)
  {
      return AccessPoint<P, D>::getDegree(p);
  }
}

namespace geolib {

using namespace PointTrait;

/*! \mainpage The GeoLib header-only library
 *
 * \section intro_sec Introduction
 *
 * This library provides functions for calculating the geographical
 * distance between two points on the Earth's surface.
 *
 * The user has the ability to define their custom point type to be
 * used for distance computation. We employ traits for accessing
 * structure elements. Therefore, these have to be specialized
 * by the library user. An example structure is shown below:
 *
 * ```cpp
 * struct CustomPoint
 * {
 *   CustomPoint(double latitude, double longitude) :
 *     latitude(latitude),
 *     longitude(longitude) {}
 *   double latitude, longitude;
 *  };
 *  ```
 *  This then has to be specialized with the generic functions
 *  `getRadian` and `getDegree`. We have to implement these in
 *  the `PointTrait` namespace. An example specialization is
 *  provided below:
 *
 * ```cpp
 * namespace PointTrait
 * {
 *     template <>
 *     struct AccessPoint<CustomPoint, 0>
 *     {
 *         static double getRadian(CustomPoint const& p)
 *         { return p.latitude * M_PI / 180; }
 *
 *         static double getDegree(CustomPoint const& p)
 *         { return p.latitude; }
 *     };
 *     template <>
 *     struct AccessPoint<CustomPoint, 1>
 *     {
 *         static double getRadian(CustomPoint const& p)
 *         { return p.longitude * M_PI / 180; }
 *
 *         static double getDegree(CustomPoint const& p)
 *         { return p.longitude; }
 *     };
 * }
 * ```
 *
 * To see this in action, please refer to the tests directory.
 *
 * \section install_sec Installation
 *
 * First clone the repository with the following command:
 * ```bash
 * $ git clone https://github.com/adl1995/boost-geometry-proposal.git
 * ```
 * Then, `cd` into the cloned repository by:
 * ```bash
 * $ cd boost-geometry-proposal
 * $ cd geolib
 * ```
 * To compile the library, first create a build directory:
 * ```bash
 * $ mkdir build
 * $ cd build
 * ```
 * The next step is to run CMake to configure the project:
 * ```bash
 * $ cmake ../
 * ```
 * Once CMake is configured, the library can be built by typing `make`. This will build the 'geolib_tests' component:
 * ```bash
 * $ make
 * ```
 * The tests can be run by typing:
 * ```bash
 * $ make tests
 * ```
 * To build the documentation using Doxygen, type:
 * ```bash
 * $ make docs
 * ```
 * Finally, the HTML documentation can be opened with Firefox by typing:
 * ```bash
 * $ firefox docs/html/index.html
 * ```
 */

/**
 * Find the geographical distance using the Haversine formula,
 * which is given by the following equation:
 *
 * \f[
 * \operatorname {hav} \left({\frac {d}{r}}\right)=\operatorname{hav}
 * (\varphi _{2}-\varphi _{1})+\cos(\varphi _{1})\cos(\varphi _{2})
 * \operatorname {hav} (\lambda _{2}-\lambda _{1})
 * \f]
 *
 * where \f${hav}\f$ is the haversine function, given by:
 *
 * \f[
 * \operatorname {hav} (\theta )=\sin ^{2}\left({\frac
 * {\theta }{2}}\right)={\frac {1-\cos(\theta )}{2}}
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Haversine_formula
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double HaversineDistance(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * Spherical Law of Cosines relates to the sides and angles of
 * spherical triangles. It is given by the following equation:
 *
 * \f[
 * d=\arccos(\sin \phi 1 * \sin \phi 2 + \cos \phi 1 *
 * \cos \phi 2 * \cos \Delta\lambda) * R
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Spherical_law_of_cosines
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double SphericalLawOfCosines(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * For less accurate and quick results, Equirectangular
 * approximation provides a good alternative. For small
 * distances Pythagoras’ theorem is used on an
 * Equi­rectangular projec­tion.
 *
 * \f[
 * {\begin{aligned}
 * x&=(\lambda -\lambda _{0})\cos \varphi _{1}\\
 * y&=(\varphi -\varphi _{1})\\
 * d&= R * \sqrt{x^{2} + y^{2}}
 * \end{aligned}}
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Equirectangular_projection
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double EquirectangularApproximation(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * The FCC (Federal Communications Commission) prescribes the
 * following formulae for distances not exceeding
 * 475 kilometres (295 mi):
 *
 * \f[
 * D={\sqrt  {(K_{1}\Delta \phi )^{2}+(K_{2}\Delta \lambda )^{2}}}
 * \f]
 *
 * where \f${K1}\f$ and \f${K2}\f$ are given by:
 *
 * \f[
 * {\begin{aligned}
 * K_{1}&=111.13209-0.56605\cos(2\phi _{m})+0.00120\cos(4\phi _{m});\\
 * K_{2}&=111.41513\cos(\phi _{m})-0.09455\cos(3\phi _{m})+
 * 0.00012\cos(5\phi _{m}).
 * \end{aligned}}
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double EllipsoidalApproximation(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * A tunnel between points on Earth is defined by a line through
 * three-dimensional space between the points of interest. This
 * is given by the following equations:
 *
 * \f[
 * {\begin{aligned}
 * &\Delta {X}=\cos(\phi _{2})\cos(\lambda _{2})-\cos(\phi _{1})\cos(\lambda _{1});\\
 * &\Delta {Y}=\cos(\phi _{2})\sin(\lambda _{2})-\cos(\phi _{1})\sin(\lambda _{1});\\
 * &\Delta {Z}=\sin(\phi _{2})-\sin(\phi _{1});\\
 * &C_{h}={\sqrt  {(\Delta {X})^{2}+(\Delta {Y})^{2}+(\Delta {Z})^{2}}}.
 * \end{aligned}}
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Geographical_distance#Tunnel_distance
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double TunnelDistance(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * Vincenty's formulae are two related iterative methods used to calculate
 * the distance between two points on the surface of a spheroid. They
 * assume an oblate spheroid model of the Earth, and hence are more
 * accurate than methods that assume a spherical Earth, such as
 * the great-circle distance.
 *
 * The inverse problem is given by:
 *
 * \f[
 * {\begin{aligned}
 * &\sin \sigma ={\sqrt  {(\cos U_{2}\sin \lambda )^{2}+(\cos U_{1}\sin U_{2}-\sin U_{1}\cos U_{2}\cos \lambda )^{2}}}\\
 * &\cos \sigma =\sin U_{1}\sin U_{2}+\cos U_{1}\cos U_{2}\cos \lambda\\
 * &\sigma =\arctan {\frac  {\sin \sigma }{\cos \sigma }}\\
 * &\sin \alpha ={\frac  {\cos U_{1}\cos U_{2}\sin \lambda }{\sin \sigma }}\\
 * &\cos ^{2}\alpha =1-\sin ^{2}\alpha\\
 * &\cos(2\sigma _{m})=\cos \sigma -{\frac  {2\sin U_{1}\sin U_{2}}{\cos ^{2}\alpha }}\\
 * &C={\frac  {f}{16}}\cos ^{2}\alpha {\big [}4+f(4-3\cos ^{2}\alpha ){\big ]}\\
 * &\lambda =L+(1-C)f\sin \alpha \left\{\sigma +C\sin \sigma \left[\cos(2\sigma _{m})+C\cos \sigma (-1+2\cos ^{2}(2\sigma _{m}))\right]\right\}\\
 * \end{aligned}}\\
 * \f]
 *
 * For more information, please refer to:
 * https://en.wikipedia.org/wiki/Vincenty's_formulae
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
double VincentysFormula(P1 const& geoPoint1, P2 const& geoPoint2);

/**
 * The distance function used in GeographicLib library. This
 * correctly handles the case for nearly antipodal points
 * and provides an accurate result upto 15 nanometers.
 *
 * For more information, please refer to:
 * https://geographiclib.sourceforge.io
 *
 * @tparam P1 The user defined point type which contains
 *     latitude and longitude.
 * @tparam P2 The user defined point type which contains
 *     latitude and longitude.
 */
template <typename P1, typename P2>
inline double GeographicLibGeodesic(P1 & geoPoint1, P2 & geoPoint2);


#endif

#ifndef GEOLIB_SRC_GEODISTANCE_HPP_IMPL
#define GEOLIB_SRC_GEODISTANCE_HPP_IMPL

template <typename P1, typename P2>
inline double HaversineDistance(P1 const& geoPoint1, P2 const& geoPoint2)
{
  // Create a local copy for convenience.
  double differenceLatitude = getRadian<0>(geoPoint2) - getRadian<0>(geoPoint1);
  double differenceLongitude = getRadian<1>(geoPoint2) - getRadian<1>(geoPoint1);

  // Calculate the central angle, given in radians.
  double centralAngle = std::pow(std::sin(differenceLatitude / 2), 2) +
      std::cos(getRadian<0>(geoPoint1)) * std::cos(getRadian<0>(geoPoint2)) *
      std::pow(std::sin(differenceLongitude / 2), 2);

  // Calculate the distance.
  double distance = (2 * std::atan2(
      std::sqrt(centralAngle), std::sqrt(1 - centralAngle)));

  return distance;
}

template <typename P1, typename P2>
inline double SphericalLawOfCosines(P1 const& geoPoint1, P2 const& geoPoint2)
{
  // Create a local copy for convenience.
  double differenceLongitude = getRadian<1>(geoPoint2) - getRadian<1>(geoPoint1);

  // Calculate the distance.
  double distance = std::acos(std::sin(getRadian<0>(geoPoint1)) *
    std::sin(getRadian<0>(geoPoint2)) + std::cos(getRadian<0>(geoPoint1)) *
    std::cos(getRadian<0>(geoPoint2)) * std::cos(differenceLongitude));

  return distance;
}

template <typename P1, typename P2>
inline double EquirectangularApproximation(P1 const& geoPoint1, P2 const& geoPoint2)
{
  // Create a local copy for convenience.
  double differenceLatitude = getRadian<0>(geoPoint2) - getRadian<0>(geoPoint1);
  double additionLatitude = getRadian<0>(geoPoint1) + getRadian<0>(geoPoint2);
  double differenceLongitude = getRadian<1>(geoPoint2) - getRadian<1>(geoPoint1);

  double x = differenceLongitude * std::cos(additionLatitude / 2);
  double y = differenceLatitude;

  // Calculate the distance.
  double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2));

  return distance;
}

template <typename P1, typename P2>
inline double EllipsoidalApproximation(P1 const& geoPoint1, P2 const& geoPoint2)
{
  // Create a local copy for convenience.
  double differenceLatitude = getDegree<0>(geoPoint2) - getDegree<0>(geoPoint1);
  double additionLatitude = getDegree<0>(geoPoint1) + getDegree<0>(geoPoint2);
  double differenceLongitude = getDegree<1>(geoPoint2) - getDegree<1>(geoPoint1);

  double K1 = 111.13209 - 0.56605 * std::cos(2 * (additionLatitude)) +
      0.00120 * std::cos(4 * (additionLatitude));
  double K2 = 111.41513 * std::cos(additionLatitude) - 0.09455 *
      std::cos(3 * (additionLatitude)) + 0.00012 *
      std::cos(5 * (additionLatitude));

  double distance = std::sqrt(std::pow(K1 * differenceLatitude, 2) +
      std::pow(K2 * differenceLongitude, 2));

  return distance;
}

template <typename P1, typename P2>
inline double TunnelDistance(P1 const& geoPoint1, P2 const& geoPoint2)
{
  double x = std::cos(getRadian<0>(geoPoint2)) * std::cos(getRadian<1>(geoPoint2)) -
      std::cos(getRadian<0>(geoPoint1)) * std::cos(getRadian<1>(geoPoint1));

  double y = std::cos(getRadian<0>(geoPoint2)) * std::sin(getRadian<1>(geoPoint2)) -
      std::cos(getRadian<0>(geoPoint1)) * std::sin(getRadian<1>(geoPoint1));

  double z = std::sin(getRadian<0>(geoPoint2)) - std::sin(getRadian<0>(geoPoint1));

  // Calculate the distance.
  double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2) +
      std::pow(z, 2));

  return distance;
}

template <typename P1, typename P2>
inline double VincentysFormula(P1 & geoPoint1, P2 & geoPoint2)
{
  // Length of semi-major axis of the ellipsoid (radius at equator).
  double a = 6378137.0; // @todo: Define a macro

  // Flattening of the ellipsoid.
  double f = 1 / 298.257223563;

  // Length of semi-minor axis of the ellipsoid (radius at the poles).
  double b = (1 - f) * a;

  double differenceLongitude = getRadian<1>(geoPoint2) - getRadian<1>(geoPoint1);

  // Reduced latitude (latitude on the auxiliary sphere).
  double U1 = std::atan((1 - f) *
      std::tan(getRadian<0>(geoPoint1)));
  double U2 = std::atan((1 - f) *
      std::tan(getRadian<0>(geoPoint2)));

  // To avoid duplication, define all the required variables.
  double cosU1, cosU2, sinU1, sinU2, cosLambda, sinLambda;
  double sinSigma, cosSigma, sigma, sinAlpha, cosSquaredAlpha,
      cos2Sigma, C;

  double lambda = differenceLongitude, lambdaPrime = 2 * M_PI;

  while (std::abs(lambda - lambdaPrime) > 1e-12)
  {
    // For convenience.
    cosU1 = std::cos(U1);
    cosU2 = std::cos(U2);
    sinU1 = std::sin(U1);
    sinU2 = std::sin(U2);
    cosLambda = std::cos(lambda);
    sinLambda = std::sin(lambda);

    // Given the coordinates, the inverse problem is applied to find
    // the azimuths α1, α2 and the ellipsoidal distance.

    sinSigma = std::sqrt(std::pow(cosU2 * sinLambda, 2) +
        std::pow(cosU1 * sinU2 - sinU1 * cosU2 * cosLambda, 2));

    cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

    sigma = std::atan2(sinSigma, cosSigma);

    sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;

    cosSquaredAlpha = 1 - std::pow(sinAlpha, 2);

    cos2Sigma = cosSigma - 2 * sinU1 * sinU2 / cosSquaredAlpha;

    C = f / 16 * cosSquaredAlpha * (4 + f * (4 - 3 * cosSquaredAlpha));

    // lambdaPrime now holds the previous lambda.
    lambdaPrime = lambda;

    lambda = differenceLongitude + (1 - C) * f * sinAlpha *
        (sigma + C * sinSigma * (cos2Sigma + C * cosSigma *
          (1 - 2 * pow(cos2Sigma, 2))));
  }

  // When lambda has converged to the desired degree of accuracy
  // (1e-12 corresponds to approximately 0.06mm),
  // evaluate the following:

  double uSquared = cosSquaredAlpha * (std::pow(a, 2) - std::pow(b, 2)) /
      std::pow(b, 2);

  double A = 1 + uSquared / 16384 * (4096 + uSquared * (-768 + uSquared *
      (320 - 175 * uSquared)));

  double B = uSquared / 1024 * (256 + uSquared * (-128 + uSquared *
      (74 - 47 * uSquared)));

  double deltaSigma = B * sinSigma * (cos2Sigma + B / 4 * (cosSigma *
      (-1 + 2 * std::pow(cos2Sigma, 2)) - B / 6 * cos2Sigma *
      (-3 + 4 * std::pow(sinSigma, 2)) * (-3 + 4 * std::pow(cos2Sigma, 2))));

  // Calculate the ellipsoidal distance.
  double distance = b * A * (sigma - deltaSigma);

  // Calculate the azimuth (angular measurement) between the given points.
  double alpha1 = std::atan2(cosU2 * sinLambda, (cosU1 * sinU2) -
      (sinU1 * cosU2 * cosLambda));
  double alpha2 = std::atan2(cosU1 * sinLambda, (-sinU1 * cosU2) +
      (cosU1 * sinU2 * cosLambda));

  return distance;
}

template <typename P1, typename P2>
inline double GeographicLibGeodesic(P1 & geoPoint1, P2 & geoPoint2)
{
  // Length of semi-major axis of the ellipsoid (radius at equator).
  double a = 6378137.0; // @todo: Define a macro

  // Flattening of the ellipsoid.
  double f = 1 / (double(298257223563LL) / 1000000000);

  // Length of semi-minor axis of the ellipsoid (radius at the poles).
  double b = (1 - f) * a;

  // Compute longitude difference.
  // double unnormalized_lon12 = std::remainder(-getRadian<1>(geoPoint1), 360) +
  //                std::remainder(getRadian<1>(geoPoint2), 360);

  // double normalized_lon12 = std::remainder(unnormalized_lon12, 360);
  // normalized_lon12 = normalized_lon12 != -180 ? normalized_lon12 : 180;
  double normalized_lon12s, normalized_lon12 = AngDiff(getRadian<1>(geoPoint1), getRadian<1>(geoPoint2), normalized_lon12s);

  // Make longitude difference positive.
  int lonsign = normalized_lon12 >= 0 ? 1 : -1;

  double y, z = 1 / 16;
  normalized_lon12 == 0 ? 0 :
    y = std::abs(normalized_lon12);
    y = y < z ? z - (z - y) : y;
    normalized_lon12 < 0 ? -y : y;

  normalized_lon12 = lonsign * AngRound(normalized_lon12);
  normalized_lon12s = AngRound((180 - normalized_lon12) - lonsign * normalized_lon12s);

  double
    lam12 = normalized_lon12 * degree(),
    slam12, clam12;
  if (normalized_lon12 > 90) {
    sincosd(normalized_lon12s, slam12, clam12);
    clam12 = -clam12;
  } else
    sincosd(normalized_lon12, slam12, clam12);

  // If really close to the equator, treat as on equator.
  // geoPoint1.latitude = AngRound(LatFix(getRadian<0>(geoPoint1)));
  // geoPoint2.latitude = AngRound(LatFix(getRadian<0>(geoPoint2)));

  // Swap points so that point with higher (abs) latitude is point 1
  // If one latitude is a nan, then it becomes lat1.
  int swapp = std::abs(geoPoint1.latitude) < std::abs(geoPoint2.latitude) ? -1 : 1;
  if (swapp < 0) {
    lonsign *= -1;
    std::swap(geoPoint1.latitude, geoPoint2.latitude);
  }
  // Make lat1 <= 0
  int latsign = geoPoint1.latitude < 0 ? 1 : -1;
  // geoPoint1.latitude *= latsign;
  // geoPoint2.latitude *= latsign;

  // Now we have
  //   0 <= lon12 <= 180,
  //   -90 <= lat1 <= 0,
  //   lat1 <= lat2 <= -lat1.

  double sbet1, cbet1, sbet2, cbet2, s12x, m12x;

  sincosd(geoPoint1.latitude, sbet1, cbet1); sbet1 *= (1 - f);
  // Ensure cbet1 = +epsilon at poles; doing the fix on beta means that sig12
  // will be <= 2*tiny for two points at the same pole.
  norm(sbet1, cbet1); cbet1 = std::max(std::sqrt(std::numeric_limits<double>::min()), cbet1);

  sincosd(geoPoint2.latitude, sbet2, cbet2); sbet2 *= (1 - f);
  // Ensure cbet2 = +epsilon at poles
  norm(sbet2, cbet2); cbet2 = std::max(std::sqrt(std::numeric_limits<double>::min()), cbet2);

  if (cbet1 < -sbet1) {
    if (cbet2 == cbet1)
      sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
  } else {
    if (abs(sbet2) == -sbet1)
      cbet2 = cbet1;
  }

  double _ep2 = (f * (2 - f)) / ((1 - f) * (1 - f));
  double
    dn1 = sqrt(1 + _ep2 * (sbet1 * sbet1)),
    dn2 = sqrt(1 + _ep2 * (sbet2 * sbet1));

  double a12, sig12;
  int nC_ = 8;
  // index zero element of this array is unused
  double Ca[nC_];

  bool meridian = geoPoint1.latitude == -90 || slam12 == 0;

  // Not a meridian.

  // somg12 > 1 marks that it needs to be calculated

  double tiny_ = std::sqrt(std::numeric_limits<double>::min());
  double calp1, calp2, salp1, salp2, M12, M21, _f1 = (1 - f), _b = a * _f1; 
  double omg12 = 0, somg12 = 2, comg12 = 0;
  if (!meridian &&
      sbet1 == 0 &&   // and sbet2 == 0
      (f <= 0 || normalized_lon12s >= f * 180)) {

    // Geodesic runs along equator
    calp1 = calp2 = 0; salp1 = salp2 = 1;
    s12x = a * lam12;
    sig12 = omg12 = lam12 / _f1;
    m12x = _b * sin(sig12);
    M12 = M21 = cos(sig12);
    a12 = normalized_lon12 / _f1;

  } else if (!meridian) {

    // Now point1 and point2 belong within a hemisphere bounded by a
    // meridian and geodesic is neither meridional or equatorial.

    // Figure a starting point for Newton's method
    double dnm;
    sig12 = InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                         lam12, slam12, clam12,
                         salp1, calp1, salp2, calp2, dnm,
                         Ca, _ep2, f, _f1);

    unsigned outmask = OUT_MASK;
    if (sig12 >= 0) {
      // Short lines (InverseStart sets salp2, calp2, dnm)
      s12x = sig12 * _b * dnm;
      m12x = sq(dnm) * _b * sin(sig12 / dnm);
      if (outmask & GEODESICSCALE)
        M12 = M21 = std::cos(sig12 / dnm);
      a12 = sig12 / degree();
      omg12 = lam12 / (_f1 * dnm);
    } else {

      // Newton's method.  This is a straightforward solution of f(alp1) =
      // lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
      // root in the interval (0, pi) and its derivative is positive at the
      // root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
      // alp1.  During the course of the iteration, a range (alp1a, alp1b) is
      // maintained which brackets the root and with each evaluation of
      // f(alp) the range is shrunk, if possible.  Newton's method is
      // restarted whenever the derivative of f is negative (because the new
      // value of alp1 is then further from the solution) or if the new
      // estimate of alp1 lies outside (0,pi); in this case, the new starting
      // guess is taken to be (alp1a + alp1b) / 2.
      //
      // initial values to suppress warnings (if loop is executed 0 times)
      static const unsigned maxit1_ = 20;
      unsigned maxit2_ = maxit1_ + std::numeric_limits<double>::digits + 10;

      double tol0_ = std::numeric_limits<double>::epsilon(), tol2_ = std::sqrt(tol0_);
      double tolb_ = tol0_ * tol2_;
      double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0;
      unsigned numit = 0;
      // Bracketing range
      double salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
      for (bool tripn = false, tripb = false;
           numit < maxit2_;
           ++numit) {
        // the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
        // WGS84 and random input: mean = 2.85, sd = 0.60
        double dv;
        double v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                          slam12, clam12,
                          salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
                          eps, domg12, numit < maxit1_, dv, Ca, _ep2, tiny_, f, _f1);
          // Reversed test to allow escape with NaNs
          if (tripb || !(std::abs(v) >= (tripn ? 8 : 1) * tol0_)) break;
          // Update bracketing values
          if (v > 0 && (numit > maxit1_ || calp1/salp1 > calp1b/salp1b))
            { salp1b = salp1; calp1b = calp1; }
          else if (v < 0 && (numit > maxit1_ || calp1/salp1 < calp1a/salp1a))
            { salp1a = salp1; calp1a = calp1; }
          if (numit < maxit1_ && dv > 0) {
            double
              dalp1 = -v/dv;
            double
              sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
              nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
            if (nsalp1 > 0 && std::abs(dalp1) < pi()) {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1;
              salp1 = nsalp1;
              norm(salp1, calp1);
              // In some regimes we don't get quadratic convergence because
              // slope -> 0.  So use convergence conditions based on epsilon
              // instead of sqrt(epsilon).
              tripn = std::abs(v) <= 16 * tol0_;
              continue;
            }
          }
          // Either dv was not positive or updated value was outside legal
          // range.  Use the midpoint of the bracket as the next estimate.
          // This mechanism is not needed for the WGS84 ellipsoid, but it does
          // catch problems with more eccentric ellipsoids.  Its efficacy is
          // such for the WGS84 test set with the starting guess set to alp1 =
          // 90deg:
          // the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
          // WGS84 and random input: mean = 4.74, sd = 0.99
          salp1 = (salp1a + salp1b)/2;
          calp1 = (calp1a + calp1b)/2;
          norm(salp1, calp1);
          tripn = false;
          tripb = (std::abs(salp1a - salp1) + (calp1a - calp1) < tolb_ ||
                   std::abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
        }
        {
          double dummy;
          // Ensure that the reduced length and geodesic scale are computed in
          // a "canonical" way, with the I2 integral.
          unsigned lengthmask = outmask |
            (outmask & (REDUCEDLENGTH | GEODESICSCALE) ? DISTANCE : NONE);
          Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                  cbet1, cbet2, lengthmask, s12x, m12x, dummy, M12, M21, Ca, nC2_, f);
        }
        m12x *= _b;
        s12x *= _b;
        a12 = sig12 / degree();
        if (outmask & AREA) {
          // omg12 = lam12 - domg12
          double sdomg12 = sin(domg12), cdomg12 = std::cos(domg12);
          somg12 = slam12 * cdomg12 - clam12 * sdomg12;
          comg12 = clam12 * cdomg12 + slam12 * sdomg12;
        }
      }
    }

    unsigned outmask = OUT_MASK;
    double s12, m12, S12, _e2 = f * (2 - f);
    double _c2 = ((sq(a) + sq(_b) *
           (_e2 == 0 ? 1 :
            eatanhe(double(1), (f < 0 ? -1 : 1) * std::sqrt(std::abs(_e2))) / _e2))
          / 2); // authalic radius squared
    int nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER;

    if (outmask & DISTANCE)
      s12 = 0 + s12x;           // Convert -0 to 0

    if (outmask & REDUCEDLENGTH)
      m12 = 0 + m12x;           // Convert -0 to 0

    if (outmask & AREA) {
      double
        // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * cbet1,
        calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0
      double alp12;
      if (calp0 != 0 && salp0 != 0) {
        double
          // From Lambda12: tan(bet) = tan(sig) * cos(alp)
          ssig1 = sbet1, csig1 = calp1 * cbet1,
          ssig2 = sbet2, csig2 = calp2 * cbet2,
          k2 = sq(calp0) * _ep2,
          eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2),
          // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
          A4 = sq(a) * calp0 * salp0 * _e2;
        norm(ssig1, csig1);
        norm(ssig2, csig2);
        C4f(eps, Ca);
        double
          B41 = SinCosSeries(false, ssig1, csig1, Ca, nC4_),
          B42 = SinCosSeries(false, ssig2, csig2, Ca, nC4_);
        S12 = A4 * (B42 - B41);
      } else
        // Avoid problems with indeterminate sig1, sig2 on equator
        S12 = 0;

      if (!meridian && somg12 > 1) {
        somg12 = sin(omg12); comg12 = cos(omg12);
      }

      if (!meridian &&
          // omg12 < 3/4 * pi
          comg12 > -double(0.7071) &&     // Long difference not too big
          sbet2 - sbet1 < double(1.75)) { // Lat difference not too big
        // Use tan(Gamma/2) = tan(omg12/2)
        // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
        // with tan(x/2) = sin(x)/(1+cos(x))
        double domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
        alp12 = 2 * std::atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                           domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
      } else {
        // alp12 = alp2 - alp1, used in atan2 so no need to normalize
        double
          salp12 = salp2 * calp1 - calp2 * salp1,
          calp12 = calp2 * calp1 + salp2 * salp1;
        // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        // salp12 = -0 and alp12 = -180.  However this depends on the sign
        // being attached to 0 correctly.  The following ensures the correct
        // behavior.
        if (salp12 == 0 && calp12 < 0) {
          salp12 = tiny_ * calp1;
          calp12 = -1;
        }
        alp12 = std::atan2(salp12, calp12);
      }
      S12 += _c2 * alp12;
      S12 *= swapp * lonsign * latsign;
      // Convert -0 to 0
      S12 += 0;
    }
    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0) {
      std::swap(salp1, salp2);
      std::swap(calp1, calp2);
      if (outmask & GEODESICSCALE)
        std::swap(M12, M21);
    }

    salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
    salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

    // Returned value in [0, 180]
    return s12;
}

#endif

} // namespace geolib
