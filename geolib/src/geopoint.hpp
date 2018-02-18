/**
 * @file geopoint.hpp
 *
 * A representation of a geographical point on the Earth's surface.
 * The coordinates are represented in latitude and longitude.
 *
 */
#ifndef GEOLIB_SRC_GEOPOINT_HPP
#define GEOLIB_SRC_GEOPOINT_HPP

#include <cmath>

namespace geolib {

/**
 * This class represents a geographical point on the Earth's surface.
 */
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
  GeoPoint(const double latitude,
           const double longitude);

  //! Get the latitude in Radians.
  double LatitudeRad() const { return latitude * M_PI / 180; }
  //! Get the longitude in Radians.
  double LongitudeRad() const { return longitude * M_PI / 180; }

  //! Get the latitude in Degrees.
  double Latitude() const { return latitude; }
  //! Modify the latitude.
  double& Latitude() { return latitude; }

  //! Get the longitude in Degrees.
  double Longitude() const { return longitude; }
  //! Modify the longitude.
  double& Longitude() { return longitude; }

 private:
  //! Geographic coordinate that specifies the north-south
  //! position of a point on the Earth's surface.
  double latitude;

  //! Geographic coordinate that specifies the east-west
  //! position of a point on the Earth's surface.
  double longitude;
};

// Define an alias to GeoPoint.
using Point = GeoPoint;

} // namespace geolib

// Include implementation.
#include "geopoint.cpp"

#endif
