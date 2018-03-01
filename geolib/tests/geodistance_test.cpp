/**
 * @file geodistance_test.cpp
 * @author Adeel Ahmad
 *
 * Test file for GeoDistance.
 *
 */

#include <boost/geometry.hpp>
#include "../src/geodistance.hpp"

#define BOOST_TEST_MODULE GeoDistanceTest
#include <boost/test/included/unit_test.hpp>

#define E_RADIUS_KM 6373
#define E_RADIUS_M 6373000

using namespace geolib;
using namespace boost::geometry;

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

  // Convert the distance to kilometers.
  double d = E_RADIUS_KM * distance.HaversineDistance();

  BOOST_TEST(d == 0.549328, boost::test_tools::tolerance(1e-3));
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

  // Convert the distance to kilometers.
  double d = E_RADIUS_KM * distance.SphericalLawOfCosines();

  BOOST_TEST(d == 0.549328, boost::test_tools::tolerance(1e-3));
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

  // Convert the distance to kilometers.
  double d = E_RADIUS_KM * distance.EquirectangularApproximation();

  BOOST_TEST(d == 0.549328, boost::test_tools::tolerance(1e-3));
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

  // Convert the distance to kilometers.
  double d = distance.EllipsoidalApproximation();

  BOOST_TEST(d == 0.522851, boost::test_tools::tolerance(1e-3));
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

  // Convert the distance to kilometers.
  double d = E_RADIUS_KM * distance.TunnelDistance();

  BOOST_TEST(d == 0.549328, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using
 * Vincentys formula.
 */
BOOST_AUTO_TEST_CASE(VincentysFormulaTest)
{
  GeoPoint point1(23.205402, 120.335066);
  GeoPoint point2(23.202188, 120.339733);

  GeoDistance distance(point1, point2);

  // The distance is returned in meters.
  double d = distance.VincentysFormula();

  BOOST_TEST(d == 595.768367, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using Boost Geometry.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryDefaultStrategyTest)
{
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;

  // Boost Geometry takes arguments in (latitude, longitude) form.
  spherical_point point1(120.335066, 23.205402);
  spherical_point point2(120.339733, 23.202188);

  // The distance is returned in meters.
  double d = E_RADIUS_M * distance(point1, point2);

  BOOST_TEST(d == 595.768367, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Thomas strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryThomasStrategyTest)
{
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;

  typedef srs::spheroid<double> stype;
  typedef strategy::distance::thomas<stype> thomas_type;

  // Boost Geometry takes arguments in (latitude, longitude) form.
  spherical_point point1(120.335066, 23.205402);
  spherical_point point2(120.339733, 23.202188);

  // The distance is returned in meters.
  double d = distance(point1, point2, thomas_type());

  BOOST_TEST(d == 595.768367, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Vincenty strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryVincentyStrategyTest)
{
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;

  typedef srs::spheroid<double> stype;
  typedef strategy::distance::vincenty<stype> vincenty_type;

  // Boost Geometry takes arguments in (latitude, longitude) form.
  spherical_point point1(120.335066, 23.205402);
  spherical_point point2(120.339733, 23.202188);

  // The distance is returned in meters.
  double d = distance(point1, point2, vincenty_type());

  BOOST_TEST(d == 595.768367, boost::test_tools::tolerance(1e-3));
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Andoyer strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryAndoyerStrategyTest)
{
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;

  typedef  srs::spheroid<double> stype;
  typedef  strategy::distance::andoyer<stype> andoyer_type;

  // Boost Geometry takes arguments in (latitude, longitude) form.
  spherical_point point1(120.335066, 23.205402);
  spherical_point point2(120.339733, 23.202188);

  // The distance is returned in meters.
  double d = distance(point1, point2, andoyer_type());

  BOOST_TEST(d == 595.768367, boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_SUITE_END();
