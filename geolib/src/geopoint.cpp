/**
 * @file geopoint.cpp
 *
 * A representation of a geographical point on the Earth's surface.
 * The coordinates are represented in latitude and longitude.
 *
 */
#ifndef GEOLIB_SRC_GEOPOINT_CPP
#define GEOLIB_SRC_GEOPOINT_CPP

#include <cmath>

// In case it hasn't been included yet.
#include "geopoint.hpp"

namespace geolib {

GeoPoint::GeoPoint(const double latitude,
                   const double longitude) :
    latitude(latitude * M_PI / 180),
    longitude(longitude * M_PI / 180)
{ /* Nothing to do. */ }

} // namespace giolib

#endif
