/**
 * @file geographiclib-bitmask.hpp
 * The code below are part of the
 * GeographicLib library
 * https://geographiclib.sourceforge.io
 */

enum captype {
  CAP_NONE = 0U,
  CAP_C1   = 1U<<0,
  CAP_C1p  = 1U<<1,
  CAP_C2   = 1U<<2,
  CAP_C3   = 1U<<3,
  CAP_C4   = 1U<<4,
  CAP_ALL  = 0x1FU,
  CAP_MASK = CAP_ALL,
  OUT_ALL  = 0x7F80U,
  OUT_MASK = 0xFF80U,       // Includes LONG_UNROLL
};


/**
* Bit masks for what calculations to do.  These masks do double duty.
* They signify to the GeodesicLine::GeodesicLine constructor and to
* Geodesic::Line what capabilities should be included in the GeodesicLine
* object.  They also specify which results to return in the general
* routines Geodesic::GenDirect and Geodesic::GenInverse routines.
* GeodesicLine::mask is a duplication of this enum.
**********************************************************************/
enum mask {
/**
 * No capabilities, no output.
 * @hideinitializer
 **********************************************************************/
NONE          = 0U,
/**
 * Calculate latitude \e lat2.  (It's not necessary to include this as a
 * capability to GeodesicLine because this is included by default.)
 * @hideinitializer
 **********************************************************************/
LATITUDE      = 1U<<7  | CAP_NONE,
/**
 * Calculate longitude \e lon2.
 * @hideinitializer
 **********************************************************************/
LONGITUDE     = 1U<<8  | CAP_C3,
/**
 * Calculate azimuths \e azi1 and \e azi2.  (It's not necessary to
 * include this as a capability to GeodesicLine because this is included
 * by default.)
 * @hideinitializer
 **********************************************************************/
AZIMUTH       = 1U<<9  | CAP_NONE,
/**
 * Calculate distance \e s12.
 * @hideinitializer
 **********************************************************************/
DISTANCE      = 1U<<10 | CAP_C1,
/**
 * Allow distance \e s12 to be used as input in the direct geodesic
 * problem.
 * @hideinitializer
 **********************************************************************/
DISTANCE_IN   = 1U<<11 | CAP_C1 | CAP_C1p,
/**
 * Calculate reduced length \e m12.
 * @hideinitializer
 **********************************************************************/
REDUCEDLENGTH = 1U<<12 | CAP_C1 | CAP_C2,
/**
 * Calculate geodesic scales \e M12 and \e M21.
 * @hideinitializer
 **********************************************************************/
GEODESICSCALE = 1U<<13 | CAP_C1 | CAP_C2,
/**
 * Calculate area \e S12.
 * @hideinitializer
 **********************************************************************/
AREA          = 1U<<14 | CAP_C4,
/**
 * Unroll \e lon2 in the direct calculation.
 * @hideinitializer
 **********************************************************************/
LONG_UNROLL   = 1U<<15,
/**
 * All capabilities, calculate everything.  (LONG_UNROLL is not
 * included in this mask.)
 * @hideinitializer
 **********************************************************************/
ALL           = OUT_ALL| CAP_ALL,
};