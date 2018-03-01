/**
 * @file geodistance.cpp
 * @author Adeel Ahmad
 *
 * A representation of a geographical point on the Earth's surface.
 * The coordinates are represented in latitude and longitude.
 *
 */
#ifndef GEOLIB_SRC_GEODISTANCE_CPP
#define GEOLIB_SRC_GEODISTANCE_CPP

// In case it hasn't been included yet.
#include "geodistance.hpp"
#include <cmath>

namespace geolib {

GeoPoint::GeoPoint(const double& latitude,
                   const double& longitude) :
    latitude(latitude),
    longitude(longitude)
{ /* Nothing to do. */ }

GeoDistance::GeoDistance(const GeoPoint& point1,
                   const GeoPoint& point2) :
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

  // Calculate the distance.
  double distance = (2 * std::atan2(
      std::sqrt(centralAngle), std::sqrt(1 - centralAngle)));

  return distance;
}

// Find the distance using Spherical Law of Cosines.
double GeoDistance::SphericalLawOfCosines()
{
  // Create a local copy for convenience.
  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  // Calculate the distance.
  double distance = std::acos(std::sin(point1.LatitudeRad()) *
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

  // Calculate the distance.
  double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2));

  return distance;
}

// Find the distance using Ellipsoidal approximation.
double GeoDistance::EllipsoidalApproximation()
{
  // Create a local copy for convenience.
  double differenceLatitude = point2.LatitudeDeg() - point1.LatitudeDeg();
  double additionLatitude = point1.LatitudeDeg() + point2.LatitudeDeg();
  double differenceLongitude = point2.LongitudeDeg() - point1.LongitudeDeg();

  double K1 = 111.13209 - 0.56605 * std::cos(2 * (additionLatitude)) +
      0.00120 * std::cos(4 * (additionLatitude));
  double K2 = 111.41513 * std::cos(additionLatitude) - 0.09455 *
      std::cos(3 * (additionLatitude)) + 0.00012 *
      std::cos(5 * (additionLatitude));

  double distance = std::sqrt(std::pow(K1 * differenceLatitude, 2) +
      std::pow(K2 * differenceLongitude, 2));

  return distance;
}

// Find the distance using tunnel distance formula.
double GeoDistance::TunnelDistance()
{
  double x = std::cos(point2.LatitudeRad()) * std::cos(point2.LongitudeRad()) -
      std::cos(point1.LatitudeRad()) * std::cos(point1.LongitudeRad());

  double y = std::cos(point2.LatitudeRad()) * std::sin(point2.LongitudeRad()) -
      std::cos(point1.LatitudeRad()) * std::sin(point1.LongitudeRad());

  double z = std::sin(point2.LatitudeRad()) - std::sin(point1.LatitudeRad());

  // Calculate the distance.
  double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2) +
      std::pow(z, 2));

  return distance;
}

// Find the distance using Vincentys formula.
double GeoDistance::VincentysFormula()
{
  // Length of semi-major axis of the ellipsoid (radius at equator).
  double a = 6378137.0; // @todo: Define a macro

  // Flattening of the ellipsoid.
  double f = 1 / 298.257223563;

  // Length of semi-minor axis of the ellipsoid (radius at the poles).
  double b = (1 - f) * a;

  double differenceLongitude = point2.LongitudeRad() - point1.LongitudeRad();

  // Reduced latitude (latitude on the auxiliary sphere).
  double U1 = std::atan((1 - f) *
      std::tan(point1.LatitudeRad()));
  double U2 = std::atan((1 - f) *
      std::tan(point2.LatitudeRad()));

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

} // namespace geolib

#endif
