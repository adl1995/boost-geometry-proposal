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
#include <numeric>
#include <boost/chrono.hpp>
#include <boost/geometry.hpp>
#include "../src/geodistance.hpp"

#include <boost/test/included/unit_test.hpp>

using namespace geolib;
using namespace boost::chrono;
using namespace boost::geometry;
using namespace boost::test_tools;

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

      GeoPoint point1user(geoData[0], geoData[1]),
               point2user(geoData[3], geoData[4]);

      distanceTestData.push_back(geoData[6]);

      globalPoints.push_back(
          std::make_pair(point1user, point2user));
    }
  }

  // Objects for file reading.
  std::string line;
  std::ifstream infile;

  // This will hold the resulting distance.
  std::vector<double> distanceTestData;

  // This will hold the execution time for each run.
  std::vector<double> time;

  // This pair holds the user define point structure used for
  // distance calculation.
  std::vector<std::pair<GeoPoint, GeoPoint>> globalPoints;

  // These object are required for Boot Geometry tests.
  typedef model::point
      <double, 2, cs::spherical_equatorial
      <degree>> spherical_point;
  typedef srs::spheroid<double> stype;

  // Used for benchmarking the execution time.
  high_resolution_clock::time_point start;
  high_resolution_clock::time_point end;
};

BOOST_FIXTURE_TEST_SUITE(GeoDistanceTest, InitTests);

/**
 * Test case for the Haversine formula.
 */
BOOST_AUTO_TEST_CASE(HaverineTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d = E_RADIUS_M * HaversineDistance(globalPoints[i].first, globalPoints[i].second);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(1e-1));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "HaverineTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case for the
 * Spherical Law of Cosines.
 */
BOOST_AUTO_TEST_CASE(SphericalLawOfCosinesTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d = E_RADIUS_M * SphericalLawOfCosines(globalPoints[i].first, globalPoints[i].second);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(1e-1));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "SphericalLawOfCosinesTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case for the
 * Equirectangular approximation.
 */
BOOST_AUTO_TEST_CASE(EquirectangularApproximationTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d = E_RADIUS_M * EquirectangularApproximation(globalPoints[i].first, globalPoints[i].second);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(0.70));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "EquirectangularApproximationTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case for the
 * Ellipsoidal approximation.
 */
BOOST_AUTO_TEST_CASE(EllipsoidalApproximationTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d =  EllipsoidalApproximation(globalPoints[i].first, globalPoints[i].second) * 1000;

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(65.6));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "EllipsoidalApproximationTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case for the
 * Tunnel distance.
 */
BOOST_AUTO_TEST_CASE(TunnelDistanceTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d = E_RADIUS_M * TunnelDistance(globalPoints[i].first, globalPoints[i].second);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(0.57));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "TunnelDistanceTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case for the
 * Vincentys formula.
 */
BOOST_AUTO_TEST_CASE(VincentysFormulaTest)
{
  for (size_t i = 0; i < distanceTestData.size(); ++i)
  {
    start = high_resolution_clock::now();

    // The distance is retuned in meters.
    double d = VincentysFormula(globalPoints[i].first, globalPoints[i].second);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == distanceTestData[i], tolerance(1e-9));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "VincentysFormulaTest: average execution time (seconds): " << avgExecTime << std::endl;
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

    start = high_resolution_clock::now();

    // Convert the distance to meters.
    double d = E_RADIUS_M * distance(point1, point2);

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == geoData[6], tolerance(1e-2));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "BoostGeometryDefaultStrategyTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case using Boost Geometry
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

    start = high_resolution_clock::now();

    // The distance is returned in meters.
    double d = distance(point1, point2, thomas_type());

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == geoData[6], tolerance(1.5e-4));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "BoostGeometryThomasStrategyTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case using Boost Geometry
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

    start = high_resolution_clock::now();

    // The distance is returned in meters.
    double d = distance(point1, point2, vincenty_type());

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == geoData[6], tolerance(1e-5));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "BoostGeometryVincentyStrategyTest: average execution time (seconds): " << avgExecTime << std::endl;
}

/**
 * Test case using Boost Geometry
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

    start = high_resolution_clock::now();

    // The distance is returned in meters.
    double d = distance(point1, point2, andoyer_type());

    end = high_resolution_clock::now();
    time.push_back(duration_cast<duration<double>>(end - start).count());

    BOOST_TEST(d == geoData[6], tolerance(1e-3));
  }

  double avgExecTime = std::accumulate(time.cbegin(), time.cend(), 0.0);
  std::cout << "BoostGeometryAndoyerStrategyTest: average execution time (seconds): " << avgExecTime << std::endl;
}

BOOST_AUTO_TEST_SUITE_END();
