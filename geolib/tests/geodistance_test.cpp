/**
 * @file geodistance_test.cpp
 *
 * Test file for GeoDistance.
 *
 */

#include "../src/geopoint.hpp"
#include "../src/geodistance.hpp"

#define BOOST_TEST_MODULE GeoDistanceTest
#include <boost/test/included/unit_test.hpp>

using namespace geolib;

BOOST_AUTO_TEST_SUITE(GeoDistanceTest);

/**
 * Simple test case for GeoDistance.
 */
BOOST_AUTO_TEST_CASE(HaverineTest)
{
  GeoPoint point1(38.898556, -77.037852);
  GeoPoint point2(38.897147, -77.043934);

  GeoDistance distance(point1, point2);

  double haversineDistance = distance.HaversineDistance();

  BOOST_TEST(haversineDistance == 0.549, boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_SUITE_END();
