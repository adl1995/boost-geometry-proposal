/**
 * @file geodistance_test.cpp
 * @author Adeel Ahmad
 *
 * Test file for GeoDistance.
 *
 */

#define E_RADIUS_KM 6373
#define E_RADIUS_M 6373000
#define BOOST_TEST_MODULE GeoDistanceTest

#include <fstream>
#include <boost/geometry.hpp>
#include "../src/geodistance.hpp"

#include <boost/test/included/unit_test.hpp>

using namespace geolib;
using namespace boost::geometry;
using namespace boost::test_tools;


// Define a custom point representation; to
// be used for distance computation.
struct GeoPointUser
{
  GeoPointUser(double latitude, double longitude) :
    latitude(latitude), longitude(longitude)
    { /* Nothing to do. */ }

  double latitude, longitude;
};

// Specialize the generic functions getRadian and getDegree
// for our GeoPointUser type. These functions are used for
// distance computation.
namespace PointTrait
{
    template <>
    struct access<GeoPointUser, 0>
    {
        static double getRadian(GeoPointUser const& p)
        { return p.latitude * M_PI / 180; }

        static double getDegree(GeoPointUser const& p)
        { return p.latitude; }
    };
    template <>
    struct access<GeoPointUser, 1>
    {
        static double getRadian(GeoPointUser const& p)
        { return p.longitude * M_PI / 180; }

        static double getDegree(GeoPointUser const& p)
        { return p.longitude; }
    };
}

struct InitTests
{
  // Read the test data in a std::vector.
  InitTests() :
      infile("GeodTest.dat")
  {
    // This will temporarily hold the data values.
    double dataField;

    while (std::getline(infile, line))
    {
      std::istringstream iss(line);

      // Push the space separated values in a vector.
      std::vector<double> geoData;
      while (iss >> dataField) { geoData.push_back(dataField); }

      GeoPoint<double> point1(geoData[0], geoData[1]),
               point2(geoData[3], geoData[4]);

      GeoDistance<double> distanceLocal(point1, point2);

      distanceGlobalTest.push_back(
          std::make_pair(distanceLocal, geoData[6]));
    }
  }

  // Objects for file reading.
  std::string line;
  std::ifstream infile;

  // The pair will hold the GeoDistance object and
  // the resulting distance.
  std::vector<std::pair<GeoDistance<double>, double>> distanceGlobalTest;

  // These object are required for Boot Geometry tests.
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;
  typedef srs::spheroid<double> stype;
};

BOOST_FIXTURE_TEST_SUITE(GeoDistanceTest, InitTests);

/**
 * Test case for GeoDistance using
 * the Haversine formula.
 */
BOOST_AUTO_TEST_CASE(HaverineTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // Convert the distance to meters.
    double d = E_RADIUS_M * distanceGlobalTest[i].first.HaversineDistance();

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(1e-1));
  }
}

/**
 * Test case for GeoDistance using
 * Spherical Law of Cosines.
 */
BOOST_AUTO_TEST_CASE(SphericalLawOfCosinesTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // Convert the distance to meters.
    double d = E_RADIUS_M * distanceGlobalTest[i].first.SphericalLawOfCosines();

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(1e-1));
  }
}

/**
 * Test case for GeoDistance using
 * Equirectangular approximation.
 */
BOOST_AUTO_TEST_CASE(EquirectangularApproximationTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // Convert the distance to meters.
    double d = E_RADIUS_M * distanceGlobalTest[i].first.EquirectangularApproximation();

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(0.55));
  }
}

/**
 * Test case for GeoDistance using
 * Ellipsoidal approximation.
 */
BOOST_AUTO_TEST_CASE(EllipsoidalApproximationTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // Convert the distance to meters.
    double d = distanceGlobalTest[i].first.EllipsoidalApproximation() * 1000;

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(7.8));
  }
}

/**
 * Test case for GeoDistance using
 * tunnel distance.
 */
BOOST_AUTO_TEST_CASE(TunnelDistanceTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // Convert the distance to meters.
    double d = E_RADIUS_M * distanceGlobalTest[i].first.TunnelDistance();

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(0.5));
  }
}

/**
 * Test case for GeoDistance using
 * Vincentys formula.
 */
BOOST_AUTO_TEST_CASE(VincentysFormulaTest)
{
  for (size_t i = 0; i < distanceGlobalTest.size(); ++i)
  {
    // The distance is retuned in meters.
    double d = distanceGlobalTest[i].first.VincentysFormula();

    BOOST_TEST(d == distanceGlobalTest[i].second, tolerance(1e-9));
  }
}

/**
 * Test case for GeoDistance using Boost Geometry.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryDefaultStrategyTest)
{
  // This will temporarily hold the data values.
  double dataField;

  // Reset the file buffer.
  infile.clear();
  infile.seekg(0, infile.beg);

  while (std::getline(infile, line))
  {
    std::istringstream iss(line);

    // Push the space separated values in a vector.
    std::vector<double> geoData;
    while (iss >> dataField) { geoData.push_back(dataField); }

    // Boost Geometry takes arguments in (latitude, longitude) form.
    spherical_point point1(geoData[1], geoData[0]),
                    point2(geoData[4], geoData[3]);

    // Convert the distance to meters.
    double d = E_RADIUS_M * distance(point1, point2);

    BOOST_TEST(d == geoData[6], tolerance(1e-2));
  }
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Thomas strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryThomasStrategyTest)
{
  // Define the strategy.
  typedef strategy::distance::thomas<stype> thomas_type;

  // This will temporarily hold the data values.
  double dataField;

  // Reset the file buffer.
  infile.clear();
  infile.seekg(0, infile.beg);

  while (std::getline(infile, line))
  {
    std::istringstream iss(line);

    // Push the space separated values in a vector.
    std::vector<double> geoData;
    while (iss >> dataField) { geoData.push_back(dataField); }

    // Boost Geometry takes arguments in (latitude, longitude) form.
    spherical_point point1(geoData[1], geoData[0]),
                    point2(geoData[4], geoData[3]);

    // The distance is returned in meters.
    double d = distance(point1, point2, thomas_type());

    BOOST_TEST(d == geoData[6], tolerance(1e-4));
  }
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Vincenty strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryVincentyStrategyTest)
{
  // Define the strategy.
  typedef strategy::distance::vincenty<stype> vincenty_type;

  // This will temporarily hold the data values.
  double dataField;

  // Reset the file buffer.
  infile.clear();
  infile.seekg(0, infile.beg);

  while (std::getline(infile, line))
  {
    std::istringstream iss(line);

    // Push the space separated values in a vector.
    std::vector<double> geoData;
    while (iss >> dataField) { geoData.push_back(dataField); }

    // Boost Geometry takes arguments in (latitude, longitude) form.
    spherical_point point1(geoData[1], geoData[0]),
                    point2(geoData[4], geoData[3]);

    // The distance is returned in meters.
    double d = distance(point1, point2, vincenty_type());

    BOOST_TEST(d == geoData[6], tolerance(1e-5));
  }
}

/**
 * Test case for GeoDistance using Boost Geometry
 * Andoyer strategy.
 */
BOOST_AUTO_TEST_CASE(BoostGeometryAndoyerStrategyTest)
{
  // Define the strategy.
  typedef strategy::distance::andoyer<stype> andoyer_type;

  // This will temporarily hold the data values.
  double dataField;

  // Reset the file buffer.
  infile.clear();
  infile.seekg(0, infile.beg);

  while (std::getline(infile, line))
  {
    std::istringstream iss(line);

    // Push the space separated values in a vector.
    std::vector<double> geoData;
    while (iss >> dataField) { geoData.push_back(dataField); }

    // Boost Geometry takes arguments in (latitude, longitude) form.
    spherical_point point1(geoData[1], geoData[0]),
                    point2(geoData[4], geoData[3]);

    // The distance is returned in meters.
    double d = distance(point1, point2, andoyer_type());

    BOOST_TEST(d == geoData[6], tolerance(1e-4));
  }
}

BOOST_AUTO_TEST_SUITE_END();
