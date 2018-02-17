/**
 * @file geodistance.cpp
 *
 * A representation of a geographical point on the Earth's surface.
 * The coordinates are represented in latitude and longitude.
 *
 */
#ifndef GEOLIB_SRC_GEODISTANCE_CPP
#define GEOLIB_SRC_GEODISTANCE_CPP

#define E_RADIUS	6373

// In case it hasn't been included yet.
#include "geodistance.hpp"
#include <cmath>

namespace geolib {

GeoDistance::GeoDistance(const GeoPoint point1,
                   const GeoPoint point2) :
    point1(point1),
    point2(point2)
{ /* Nothing to do. */ }

// Find the distance.
double GeoDistance::HaversineDistance()
{
  // Create a local copy for convenience.
  double differenceLatitude = point2.Latitude() - point1.Latitude();
  double differenceLongitude = point2.Longitude() - point1.Longitude();

  // Calculate the central angle, given in radians.
  double centralAngle = std::pow(std::sin(differenceLatitude / 2), 2) +
      std::cos(point1.Latitude()) * std::cos(point2.Latitude()) *
      std::pow(std::sin(differenceLongitude / 2), 2);

  double distance = E_RADIUS * (2 * std::atan2(
      std::sqrt(centralAngle), std::sqrt(1 - centralAngle)));

  return distance;
}

} // namespace giolib

#endif
