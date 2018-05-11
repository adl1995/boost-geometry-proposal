#include <iostream>
#include <numeric>
#include <boost/geometry.hpp>
#include "../src/geodistance.hpp"

using namespace geolib;
using namespace boost::geometry;

/**
 * This user defined structure represents a 2-dimensional
 * geographic / spherical point on the Earth's surface
 * using degrees as units.
 */
struct GeoPoint
{
  /**
   * Construct the GeoPoint object with the given parameters.
   *
   * @param latitude  Geographic coordinate that specifies the
   *     north-south position of a point on the Earth's surface.
   * @param longitude  Geographic coordinate that specifies the
   *     east-west position of a point on the Earth's surface.
   */
  GeoPoint(double latitude, double longitude) :
    latitude(latitude),
    longitude(longitude)
  { /* Nothing to do. */ }

  void setLatitude(double lat) { latitude = lat; }
  void setLongitude(double lon) { longitude = lon; }

  public:
    double latitude, longitude;
};

/**
 * Specialize the generic functions getRadian and getDegree
 * for our GeoPoint type. These functions are used for
 * distance computation.
 */
namespace PointTrait
{
  template <>
  struct AccessPoint<GeoPoint, 0>
  {
    static double getRadian(GeoPoint const& p)
    { return p.latitude * M_PI / 180; }

    static double getDegree(GeoPoint const& p)
    { return p.latitude; }

    // TODO: add methods for updating structure members.
  };
  template <>
  struct AccessPoint<GeoPoint, 1>
  {
    static double getRadian(GeoPoint const& p)
    { return p.longitude * M_PI / 180; }

    static double getDegree(GeoPoint const& p)
    { return p.longitude; }
  };
}


/**
 * Test case for the
 * GeographicLibGeodesic.
 */
int main()
{
  GeoPoint point1user(36.530042355041, 0),
           point2user(-48.164270779097768864, 5.762344694676510456);
  // Distance (returned): 9381056.7211346551776 meters
  // Distance (actual): 9398502.0434630513191 meters

  double d = GeographicLibGeodesic(point1user, point2user);

  std::cout << std::setprecision(20);
  std::cout << "GeographicLibGeodesic: distance: " << d << std::endl;

  point1user.setLatitude(2.179167), point1user.setLongitude(-73.787500);
  point2user.setLatitude(-2.1622), point2user.setLongitude(106.139064);
  // Distance (returned): 19992302.326442237943
  // Distance (actual): 20001571

  // The distance is retuned in meters.
  d = GeographicLibGeodesic(point1user, point2user);

  std::cout << std::setprecision(20);
  std::cout << "GeographicLibGeodesic: distance: " << d << std::endl;
}