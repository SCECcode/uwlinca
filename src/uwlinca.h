/**
 * @file uwlinca.h
 * @brief Main header file for uwlinca library.
 * @version 1.0
 *
 * Delivers the uwlinca model 
 * base on original uwlinca
 *
 */

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "proj.h"

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD M_PI / 180.0

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1

/* config string */
#define UWLINCA_CONFIG_MAX 1000

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct uwlinca_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} uwlinca_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct uwlinca_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
	/** Qp */
	double qp;
	/** Qs */
	double qs;
} uwlinca_properties_t;

/** The CVM-S5 configuration structure. */
typedef struct uwlinca_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting in UTM projection */
	double top_left_corner_e;
	/** Top left corner northing in UTM projection */
	double top_left_corner_n;
	/** Top right corner easting in UTM projection */
	double top_right_corner_e;
	/** Top right corner northing in UTM projection */
	double top_right_corner_n;
	/** Bottom left corner easting in UTM projection */
	double bottom_left_corner_e;
	/** Bottom left corner northing in UTM projection */
	double bottom_left_corner_n;
	/** Bottom right corner easting in UTM projection */
	double bottom_right_corner_e;
	/** Bottom right corner northing in UTM projection */
	double bottom_right_corner_n;
	/** Z interval for the data */
	double depth_interval;
        /** The data access seek method, fast-X, or fast-Y */
        char seek_axis[128];
        /** The data seek direction, bottom-up, or top-down */
        char seek_direction[128];
        /** trilinear interploation; */
        int interpolation;
} uwlinca_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct uwlinca_model_t {
	/** A pointer to the Vs data either in memory or disk. Null if does not exist. */
	void *vs;
	/** Vs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vs_status;
	/** A pointer to the Vp data either in memory or disk. Null if does not exist. */
	void *vp;
	/** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vp_status;
	/** A pointer to the rho data either in memory or disk. Null if does not exist. */
	void *rho;
	/** Rho status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int rho_status;
	/** A pointer to the Qp data either in memory or disk. Null if does not exist. */
	void *qp;
	/** Qp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qp_status;
	/** A pointer to the Qs data either in memory or disk. Null if does not exist. */
	void *qs;
	/** Qs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qs_status;
} uwlinca_model_t;

// UCVM API Required Functions
#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(uwlinca_point_t *points, uwlinca_properties_t *data, int numpts);
int model_config(char **config, int *sz);

int (*get_model_init())(const char *, const char *);
int (*get_model_query())(uwlinca_point_t *, uwlinca_properties_t *, int);
int (*get_model_finalize())();
int (*get_model_version())(char *, int);
int (*get_model_config())(char **, int*);

#endif

// uwlinca Related Functions

/** Initializes the model */
int uwlinca_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int uwlinca_finalize();
/** Returns version information */
int uwlinca_version(char *ver, int len);
/** Queries the model */
int uwlinca_query(uwlinca_point_t *points, uwlinca_properties_t *data, int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int uwlinca_read_configuration(char *file, uwlinca_configuration_t *config);
int uwlinca_dump_configuration(uwlinca_configuration_t *config);
/** Prints out the error string. */
void uwlinca_print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void uwlinca_read_properties(int x, int y, int z, uwlinca_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int uwlinca_try_reading_model(uwlinca_model_t *model);

// Interpolation Functions
/** Linearly interpolates two uwlinca_properties_t structures */
void uwlinca_linear_interpolation(double percent, uwlinca_properties_t *x0, uwlinca_properties_t *x1, uwlinca_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void uwlinca_bilinear_interpolation(double x_percent, double y_percent, uwlinca_properties_t *four_points, uwlinca_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void uwlinca_trilinear_interpolation(double x_percent, double y_percent, double z_percent, uwlinca_properties_t *eight_points,
							 uwlinca_properties_t *ret_properties);
