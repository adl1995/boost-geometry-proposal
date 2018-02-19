/**
 * @file geodistance_test.cpp
 *
 * Test file for GeoDistance.
 *
 */

#define BOOST_TEST_MODULE GeoDistanceTest
#include <boost/test/included/unit_test.hpp>

#include "../src/geopoint.hpp"
#include "../src/geodistance.hpp"

using namespace geolib;

BOOST_AUTO_TEST_SUITE(GeoDistanceTest);

/**
 * Test case for GeoDistance using
 * the Haversine formula.
 */
BOOST_AUTO_TEST_CASE(HaverineTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double d = distance.HaversineDistance();

  BOOST_TEST(d == 0.549, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using
 * Spherical Law of Cosines.
 */
BOOST_AUTO_TEST_CASE(SphericalLawOfCosinesTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double d = distance.SphericalLawOfCosines();

  BOOST_TEST(d == 0.549, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using
 * Equirectangular approximation.
 */
BOOST_AUTO_TEST_CASE(EquirectangularApproximationTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double d = distance.EquirectangularApproximation();

  BOOST_TEST(d == 0.549, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using
 * Ellipsoidal approximation.
 */
BOOST_AUTO_TEST_CASE(EllipsoidalApproximationTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double d = distance.EllipsoidalApproximation();

  BOOST_TEST(d == 0.549, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using
 * tunnel distance.
 */
BOOST_AUTO_TEST_CASE(TunnelDistanceTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double d = distance.TunnelDistance();

  BOOST_TEST(d == 0.549, boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_SUITE_END();
