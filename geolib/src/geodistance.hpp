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

namespace geolib {

/**
 * This class represents a 2-dimensional geographic / spherical point
 * on the Earth's surface using degrees as units.
 */
template <typename CoordType>
class GeoPoint
{
 public:
  /**
   * Construct the GeoPoint object with the given parameters.
   *
   * @param latitude  Geographic coordinate that specifies the
   *     north-south position of a point on the Earth's surface.
   * @param longitude  Geographic coordinate that specifies the
   *     east-west position of a point on the Earth's surface.
   */
  GeoPoint();
  GeoPoint(const CoordType& latitude,
           const CoordType& longitude);

  //! Get the latitude in radians.
  CoordType LatitudeRad() const { return latitude * M_PI / 180; }
  //! Get the longitude in radians.
  CoordType LongitudeRad() const { return longitude * M_PI / 180; }

  //! Get the latitude in degrees.
  CoordType LatitudeDeg() const { return latitude; }
  //! Modify the latitude.
  CoordType& Latitude() { return latitude; }

  //! Get the longitude in degrees.
  CoordType LongitudeDeg() const { return longitude; }
  //! Modify the longitude.
  CoordType& Longitude() { return longitude; }

 private:
  //! Geographic coordinate that specifies the north-south
  //! position of a point on the Earth's surface.
  CoordType latitude;

  //! Geographic coordinate that specifies the east-west
  //! position of a point on the Earth's surface.
  CoordType longitude;
};

/**
 * This class calculates the geographical distance between
 * two points on the Earth's surface.
 */
template <typename CoordType>
class GeoDistance
{
 public:
  /**
   * Construct the GeoDistance object with the given parameters.
   *
   * @param point1 GeoPoint object for the first position
   *     on the Earth's surface.
   * @param point2 GeoPoint object for the second position
   *     on the Earth's surface.
   */
  GeoDistance(const GeoPoint<CoordType>& point1,
              const GeoPoint<CoordType>& point2);

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
   */
  double HaversineDistance();

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
   */
  double SphericalLawOfCosines();

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
   */
  double EquirectangularApproximation();

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
   */
  double EllipsoidalApproximation();

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
   */
  double TunnelDistance();

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
   */
  double VincentysFormula();

  //! Get the first point.
  GeoPoint<CoordType> Point1() const { return point1; }
  //! Modify the first point.
  GeoPoint<CoordType>& Point1() { return point1; }

  //! Get the second point.
  GeoPoint<CoordType> Point2() const { return point2; }
  //! Modify the second point.
  GeoPoint<CoordType>& Point2() { return point2; }

 private:
  //! GeoPoint object for the first position on the Earth's surface.
  GeoPoint<CoordType> point1;

  //! GeoPoint object for the second position on the Earth's surface.
  GeoPoint<CoordType> point2;
};

} // namespace geolib

// Include implementation.
#include "geodistance.cpp"

#endif
