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
inline double VincentysFormula(P1 const& geoPoint1, P2 const& geoPoint2)
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

#endif

} // namespace geolib
