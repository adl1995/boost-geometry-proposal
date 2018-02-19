/**
 * @file geodistance.hpp
 *
 * Geographical distance between two points on the Earth's surface.
 * The current implementation includes the Haversine formula.
 *
 */
#ifndef GEOLIB_SRC_GEODISTANCE_HPP
#define GEOLIB_SRC_GEODISTANCE_HPP

#include "geopoint.hpp"

namespace geolib {

/**
 * This class calculates the geographical distance between
 * two points on the Earth's surface.
 */
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
  GeoDistance(const GeoPoint point1,
              const GeoPoint point2);

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
   * d = \arccos(\sin \phi 1 * \sin \phi 2 + \cos \phi 1 *
   * \cos \phi 2 * \cos \Delta\lambda) * R
   * \f]
   *
   * For more information, please refer to:
   * https://en.wikipedia.org/wiki/Spherical_law_of_cosines
   */
  double SphericalLawOfCosines();

  //! Get the first point.
  GeoPoint Point1() const { return point1; }
  //! Modify the first point.
  GeoPoint& Point1() { return point1; }

  //! Get the second point.
  GeoPoint Point2() const { return point2; }
  //! Modify the second point.
  GeoPoint& Point2() { return point2; }

 private:
  //! GeoPoint object for the first position on the Earth's surface.
  GeoPoint point1;

  //! GeoPoint object for the second position on the Earth's surface.
  GeoPoint point2;
};

} // namespace geolib

// Include implementation.
#include "geodistance.cpp"

#endif
