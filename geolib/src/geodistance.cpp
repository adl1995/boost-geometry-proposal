/**
 * @file geodistance.cpp
 *
 * A representation of a geographical point on the Earth's surface.
 * The coordinates are represented in latitude and longitude.
 *
 */
#ifndef GEOLIB_SRC_GEODISTANCE_CPP
#define GEOLIB_SRC_GEODISTANCE_CPP

#define E_RADIUS_KM	6373

// In case it hasn't been included yet.
#include "geodistance.hpp"
#include <cmath>

namespace geolib {

GeoDistance::GeoDistance(const GeoPoint point1,
                   const GeoPoint point2) :
    point1(point1),
    point2(point2)
{ /* Nothing to do. */ }

// Find the distance using the Haversine formula.
double GeoDistance::HaversineDistance()
{
  // Create a local copy for convenience.
  double differenceLatitude = point2.LatitudeRad() - point1.LatitudeRad();
  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  // Calculate the central angle, given in radians.
  double centralAngle = std::pow(std::sin(differenceLatitude / 2), 2) +
      std::cos(point1.LatitudeRad()) * std::cos(point2.LatitudeRad()) *
      std::pow(std::sin(differenceLongitude / 2), 2);

  // Calculate the distance in kilometers.
  double distance = E_RADIUS_KM * (2 * std::atan2(
      std::sqrt(centralAngle), std::sqrt(1 - centralAngle)));

  return distance;
}

// Find the distance using Spherical Law of Cosines.
double GeoDistance::SphericalLawOfCosines()
{
  // Create a local copy for convenience.
  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  // Calculate the distance in kilometers.
  double distance = E_RADIUS_KM * std::acos(std::sin(point1.LatitudeRad()) *
    std::sin(point2.LatitudeRad()) + std::cos(point1.LatitudeRad()) *
    std::cos(point2.LatitudeRad()) * std::cos(differenceLongitude));

  return distance;
}

// Find the distance using Equirectangular approximation.
double GeoDistance::EquirectangularApproximation()
{
  // Create a local copy for convenience.
  double differenceLatitude = point2.LatitudeRad() - point1.LatitudeRad();
  double additionLatitude = point1.LatitudeRad() + point2.LatitudeRad();
  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  double x = differenceLongitude * std::cos(additionLatitude / 2);
  double y = differenceLatitude;

  // Calculate the distance in kilometers.
  double distance = E_RADIUS_KM * std::sqrt(std::pow(x, 2) + std::pow(y, 2));

  return distance;
}

// Find the distance using Ellipsoidal approximation.
double GeoDistance::EllipsoidalApproximation()
{
  // Create a local copy for convenience.
  double differenceLatitude = point2.LatitudeRad() - point1.LatitudeRad();
  double additionLatitude = point1.LatitudeRad() + point2.LatitudeRad();
  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  double K1 = 111.13209 - 0.56605 * std::cos(2 * (additionLatitude / 2)) +
      0.00120 * std::cos(4 * (additionLatitude / 2));
  double K2 = 111.41513 * std::cos(additionLatitude / 2) - 0.09455 *
      std::cos(3 * (additionLatitude / 2)) + 0.00012 *
      std::cos(5 * (additionLatitude / 2));

  double distance = std::sqrt(std::pow(K1 * differenceLatitude, 2) +
      std::pow(K2 * differenceLongitude, 2));

  return distance;
}

// Find the distance using theh tunnel distance formula.
double GeoDistance::TunnelDistance()
{
  double x = std::cos(point2.LatitudeRad()) * std::cos(point2.LongitudeRad()) -
      std::cos(point1.LatitudeRad()) * std::cos(point1.LongitudeRad());

  double y = std::cos(point2.LatitudeRad()) * std::sin(point2.LongitudeRad()) -
      std::cos(point1.LatitudeRad()) * std::sin(point1.LongitudeRad());

  double z = std::sin(point2.LatitudeRad()) - std::sin(point1.LatitudeRad());

  // Calculate the distance in kilometers.
  double distance = E_RADIUS_KM * std::sqrt(std::pow(x, 2) + std::pow(y, 2) +
      std::pow(z, 2));

  return distance;
}

} // namespace geolib

#endif
