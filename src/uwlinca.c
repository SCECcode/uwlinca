/**
 * @file uwlinca.c
 * @brief Main file for uwlinca model
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the uwlinca model which base from 
 * uwlinca that was created to work with mesh created by ucvm2mesh
 *
 */

#include "limits.h"
#include "ucvm_model_dtypes.h"
#include "uwlinca.h"
#include <assert.h>

int uwlinca_debug=0;

/** The config of the model */
char *uwlinca_config_string=NULL;
int uwlinca_config_sz=0;

// Constants
/** The version of the model. */
const char *uwlinca_version_string = "uwlinca";

// Variables
/** Set to 1 when the model is ready for query. */
int uwlinca_is_initialized = 0;

char uwlinca_data_directory[128];

/** Configuration parameters. */
uwlinca_configuration_t *uwlinca_configuration;
/** Holds pointers to the velocity model data OR indicates it can be read from file. */
uwlinca_model_t *uwlinca_velocity_model;

/** Proj coordinate transformation objects. can go from geo <-> utm */
PJ *uwlinca_geo2utm = NULL;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double uwlinca_cos_rotation_angle = 0;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double uwlinca_sin_rotation_angle = 0;

/** The height of this model's region, in meters. */
double uwlinca_total_height_m = 0;
/** The width of this model's region, in meters. */
double uwlinca_total_width_m = 0;
/**
 * Initializes the uwlinca plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int uwlinca_init(const char *dir, const char *label) {
    int tempVal = 0;
    char configbuf[512];
    double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

    // Initialize variables.
    uwlinca_configuration = calloc(1, sizeof(uwlinca_configuration_t));
    uwlinca_velocity_model = calloc(1, sizeof(uwlinca_model_t));

    uwlinca_config_string = calloc(UWLINCA_CONFIG_MAX, sizeof(char));
    uwlinca_config_string[0]='\0';
    uwlinca_config_sz=0;

    // Configuration file location.
    sprintf(configbuf, "%s/model/%s/data/config", dir, label);

    // Read the uwlinca_configuration file.
    if (uwlinca_read_configuration(configbuf, uwlinca_configuration) != SUCCESS)
        return FAIL;

    // Set up the data directory.
    sprintf(uwlinca_data_directory, "%s/model/%s/data/%s/", dir, label, uwlinca_configuration->model_dir);

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    tempVal = uwlinca_try_reading_model(uwlinca_velocity_model);

    if (tempVal == SUCCESS) {
//        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
//        fprintf(stderr, "hard disk may result in slow performance.\n");
    } else if (tempVal == FAIL) {
        uwlinca_print_error("No model file was found to read from.");
        return FAIL;
    }

    // We need to convert the point from lat, lon to UTM, let's set it up.
    char uwlinca_projstr[64];
    snprintf(uwlinca_projstr, 64, "+proj=utm +ellps=clrk66 +zone=%d +datum=NAD27 +units=m +no_defs", uwlinca_configuration->utm_zone);
    if (!(uwlinca_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", uwlinca_projstr, NULL))) {
        uwlinca_print_error("Could not set up Proj transformation from EPSG:4326 to UTM.");
        uwlinca_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }


    // In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
    // corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
    // point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
    // the X and Y axis determines which grid points we use for the interpolation routine.

    // Calculate the rotation angle of the box.
    north_height_m = uwlinca_configuration->top_left_corner_n - uwlinca_configuration->bottom_left_corner_n;
    east_width_m = uwlinca_configuration->top_left_corner_e - uwlinca_configuration->bottom_left_corner_e;

    // Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
    rotation_angle = atan(east_width_m / north_height_m);

    uwlinca_cos_rotation_angle = cos(rotation_angle);
    uwlinca_sin_rotation_angle = sin(rotation_angle);

    uwlinca_total_height_m = sqrt(pow(uwlinca_configuration->top_left_corner_n - uwlinca_configuration->bottom_left_corner_n, 2.0f) +
          pow(uwlinca_configuration->top_left_corner_e - uwlinca_configuration->bottom_left_corner_e, 2.0f));
    uwlinca_total_width_m  = sqrt(pow(uwlinca_configuration->top_right_corner_n - uwlinca_configuration->top_left_corner_n, 2.0f) +
          pow(uwlinca_configuration->top_right_corner_e - uwlinca_configuration->top_left_corner_e, 2.0f));

    if(uwlinca_debug) {
      fprintf(stderr,"north_height %lf east_width %lf\n", north_height_m, east_width_m);
      fprintf(stderr,"totol height %lf total width %lf\n", uwlinca_total_height_m, uwlinca_total_width_m);
      fprintf(stderr,"cos angle %lf sin angle %lf\n", uwlinca_cos_rotation_angle, uwlinca_sin_rotation_angle);
    }

    // setup config_string 
    sprintf(uwlinca_config_string,"config = %s\n",configbuf);
    uwlinca_config_sz=1;


    // Let everyone know that we are initialized and ready for business.
    uwlinca_is_initialized = 1;

    return SUCCESS;
}

/*
#define PROJ_GEO_IPJ "+proj=latlong +datum=WGS84"
#define PROJ_GEO_OPJ "+proj=utm +zone=11 +ellps=WGS84"
static int to_utm(double *lon, double *lat) {
    projPJ ipj= pj_init_plus(PROJ_GEO_IPJ);
    projPJ opj = pj_init_plus(PROJ_GEO_OPJ);
    int p = pj_transform(ipj, opj, 1, 1, lon, lat, NULL );
    return p;
}
*/

/*** transform to UTM zone 32, then back to geographical 
         https://proj.org/en/9.3/development/quickstart.html
14    P = proj_create_crs_to_crs(
15        C, "EPSG:4326", "+proj=utm +zone=32 +datum=WGS84", 
16        NULL);

41    b = proj_trans(P, PJ_FWD, a);
42    printf("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);
43
44    b = proj_trans(P, PJ_INV, b);
45    printf("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);

and
    PJ_COORD c_in;
33    c_in.lpzt.z = 0.0;
34    c_in.lpzt.t = HUGE_VAL; // important only for time-dependent projections
35    c_in.lp.lam = lon;
36    c_in.lp.phi = lat;
37

and

typedef union {
    double v[4];
    PJ_XYZT xyzt;
    PJ_UVWT uvwt;
    PJ_LPZT lpzt;
    PJ_GEOD geod;
    PJ_OPK opk;
    PJ_ENU enu;
    PJ_XYZ  xyz;
    PJ_UVW  uvw;
    PJ_LPZ  lpz;
    PJ_XY   xy;
    PJ_UV   uv;
    PJ_LP   lp;
} PJ_COORD ;

examples:
https://proj.org/en/5.0/development/migration.html

*/

static int to_utm(double lon, double lat, double *point_u, double *point_v) {
    PJ_COORD xyzSrc = proj_coord(lat, lon, 0.0, HUGE_VAL);
    PJ_COORD xyzDest = proj_trans(uwlinca_geo2utm, PJ_FWD, xyzSrc);
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n",
              lat, lon);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *point_u = xyzDest.xyzt.x;
    *point_v = xyzDest.xyzt.y;
    return err;
}

static int to_geo(double point_u, double point_v, double *lon, double *lat) {
    PJ_COORD xyzSrc;
    xyzSrc.xyzt.x=point_u;
    xyzSrc.xyzt.y=point_v;
    PJ_COORD xyzDest = proj_trans(uwlinca_geo2utm, PJ_INV, xyzSrc);
    
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming u=%.4f, v=%.4f to Geo.\n",
              point_u, point_v);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *lon=xyzDest.lp.lam;
    *lat=xyzDest.lp.phi;
    return err;
}


/**
 * Queries uwlinca at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, rho, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int uwlinca_query(uwlinca_point_t *points, uwlinca_properties_t *data, int numpoints) {
    int i = 0;

    double point_u = 0, point_v = 0;
    double point_x = 0, point_y = 0; 
				   
    int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
    double x_percent = 0, y_percent = 0, z_percent = 0;

    uwlinca_properties_t surrounding_points[8];
    int zone = uwlinca_configuration->utm_zone;

    for (i = 0; i < numpoints; i++) {

        // We need to be below the surface to service this query.
        if (points[i].depth < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

	// lon,lat,u,v			     
	to_utm(points[i].longitude, points[i].latitude, &point_u, &point_v);

if(uwlinca_debug) { fprintf(stderr,"lon %lf lat %lf\n", points[i].longitude, points[i].latitude); }
if(uwlinca_debug) { fprintf(stderr,"point_u %lf point_v %lf\n", point_u, point_v); }

        // Point within rectangle.
        point_u -= uwlinca_configuration->bottom_left_corner_e;
        point_v -= uwlinca_configuration->bottom_left_corner_n;

        // We need to rotate that point, the number of degrees we calculated above.
        point_x = uwlinca_cos_rotation_angle * point_u - uwlinca_sin_rotation_angle * point_v;
        point_y = uwlinca_sin_rotation_angle * point_u + uwlinca_cos_rotation_angle * point_v;

        // Which point base point does that correspond to?
        load_x_coord = floor(point_x / uwlinca_total_width_m * (uwlinca_configuration->nx - 1));

/* make origin-y at lower left instead of upper left (flipped) */
        load_y_coord = floor(point_y / uwlinca_total_height_m * (uwlinca_configuration->ny - 1));
if(uwlinca_debug) { fprintf(stderr,"  before load_y_coord %d\n", load_y_coord); }
        load_y_coord = (uwlinca_configuration->ny - load_y_coord) - 1;
if(uwlinca_debug) { fprintf(stderr,"  after load_y_coord (%d)%d\n", uwlinca_configuration->ny,load_y_coord); }

        // And on the Z-axis?
        load_z_coord = (uwlinca_configuration->depth / uwlinca_configuration->depth_interval - 1) -
                       floor(points[i].depth / uwlinca_configuration->depth_interval);

if(uwlinca_debug) { fprintf(stderr,"load_x_coord %d load_y_coord %d load_z_coord %d\n", load_x_coord,load_y_coord,load_z_coord); }

        // Are we outside the model's X and Y boundaries?
        if (load_x_coord > uwlinca_configuration->nx - 2 || load_y_coord > uwlinca_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        if(uwlinca_configuration->interpolation) {

          // Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
          double x_interval=(uwlinca_configuration->nx > 1) ?
                     uwlinca_total_width_m / (uwlinca_configuration->nx-1):uwlinca_total_width_m;
          double y_interval=(uwlinca_configuration->ny > 1) ?
                     uwlinca_total_height_m / (uwlinca_configuration->ny-1):uwlinca_total_height_m;

          x_percent = fmod(point_u, x_interval) / x_interval;
          y_percent = fmod(point_v, y_interval) / y_interval;
          z_percent = fmod(points[i].depth, uwlinca_configuration->depth_interval) / uwlinca_configuration->depth_interval;

          if (load_z_coord < 1) {
              // We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
              data[i].vp = -1;
              data[i].vs = -1;
              data[i].rho = -1;
              data[i].qp = -1;
              data[i].qs = -1;
              continue;
          } else {
              // Read all the surrounding point properties.
              uwlinca_read_properties(load_x_coord, load_y_coord, load_z_coord, &(surrounding_points[0]));    // Orgin.
              uwlinca_read_properties(load_x_coord + 1, load_y_coord, load_z_coord, &(surrounding_points[1]));    // Orgin + 1x
              uwlinca_read_properties(load_x_coord, load_y_coord + 1, load_z_coord, &(surrounding_points[2]));    // Orgin + 1y
              uwlinca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord, &(surrounding_points[3]));    // Orgin + x + y, forms top plane.
              uwlinca_read_properties(load_x_coord, load_y_coord, load_z_coord - 1, &(surrounding_points[4]));    // Bottom plane origin
              uwlinca_read_properties(load_x_coord + 1, load_y_coord, load_z_coord - 1, &(surrounding_points[5]));    // +1x
              uwlinca_read_properties(load_x_coord, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));    // +1y
              uwlinca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));    // +x +y, forms bottom plane.
  
              uwlinca_trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
          }
          } else {
if(uwlinca_debug) {fprintf(stderr,"direct call, no interpolation\n"); }
              uwlinca_read_properties(load_x_coord, load_y_coord, load_z_coord, &(data[i]));    // Orgin.
        }

    }

    return SUCCESS;
}

/* Density derived from Vp via Nafe-Drake curve, Brocher (2005) eqn 1. */
double _get_rho(double f) {
  double rho;

  /* Convert m to km */
  f = f / 1000.0;
  rho = f * (1.6612 - f * (0.4721 - f * (0.0671 - f * (0.0043 - f * 0.000106))));
  if (rho < 1.0) {
    rho = 1.0;
  }
  rho = rho * 1000.0;
  return(rho);
}

**
 * Calculates the vs based off of Vp. Base on Brocher's formulae
 *
 * https://pubs.usgs.gov/of/2005/1317/of2005-1317.pdf
 *
 * @param vp
 * @return Vs, in km.
 * Vs derived from Vp, Brocher (2005) eqn 1.
 * [eqn. 1] Vs (km/s) = 0.7858 – 1.2344Vp + 0.7949Vp2 – 0.1238Vp3 + 0.0064Vp4.
 * Equation 1 is valid for 1.5 < Vp < 8 km/s.
 */
double uwlinca_calculate_vs(double vp) {
     double retVal ;

     vp = vp * 0.001;
     double t1= (vp * 1.2344);
     double t2= ((vp * vp)* 0.7949);
     double t3= ((vp * vp * vp) * 0.1238);
     double t4= ((vp * vp * vp * vp) * 0.0064);
     retVal = 0.7858 - t1 + t2 - t3 + t4;
     retVal = retVal * 1000.0;

     return retVal;
}

/**
 * Calculates the density based off of Vp. Base on Brocher's formulae
 *
 * @param vp
 * @return Density, in g/m^3.
 * [eqn. 6] r (g/cm3) = 1.6612Vp – 0.4721Vp2 + 0.0671Vp3 – 0.0043Vp4 + 0.000106Vp5.
 * Equation 6 is the “Nafe-Drake curve” (Ludwig et al., 1970).
 * start with vp in km
 */
double uwlinca_calculate_density(double vp) {
     double retVal ;

     vp = vp * 0.001;
     double t1 = (vp * 1.6612);
     double t2 = ((vp * vp ) * 0.4721);
     double t3 = ((vp * vp * vp) * 0.0671);
     double t4 = ((vp * vp * vp * vp) * 0.0043);
     double t5 = ((vp * vp * vp * vp * vp) * 0.000106);
     retVal = t1 - t2 + t3 - t4 + t5;
     if (retVal < 1.0) {
       retVal = 1.0;
     }
     retVal = retVal * 1000.0;
     return retVal;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void uwlinca_read_properties(int x, int y, int z, uwlinca_properties_t *data) {
    // Set everything to -1 to indicate not found.
    data->vp = -1;
    data->vs = -1;
    data->rho = -1;
    data->qp = -1;
    data->qs = -1;

if(uwlinca_debug) {fprintf(stderr,"read_properties index: x(%d) y(%d) z(%d)\n",x,y,z); }
if(uwlinca_debug) {fprintf(stderr,"     nx(%d) ny(%d) nz(%d)\n",
	          uwlinca_configuration->nx,uwlinca_configuration->ny,uwlinca_configuration->nz); }
    float *ptr = NULL;
    FILE *fp = NULL;
    long location = 0;

    // the z is inverted at line #145
    if ( strcmp(uwlinca_configuration->seek_axis, "fast-y") == 0 ||
                 strcmp(uwlinca_configuration->seek_axis, "fast-Y") == 0 ) { // fast-y,  uwlinca 
        if(strcmp(uwlinca_configuration->seek_direction, "bottom-up") == 0) { 
                location = ((long) z * uwlinca_configuration->nx * uwlinca_configuration->ny) + (x * uwlinca_configuration->ny) + y;
if(uwlinca_debug) {fprintf(stderr,"LOCATION==%d(fast-y, bottom-up)\n", location); }
            } else { // nz starts from 0 up to nz-1
                    location = ((long)((uwlinca_configuration->nz -1) - z) * uwlinca_configuration->nx * uwlinca_configuration->ny) + (x * uwlinca_configuration->ny) + y;
if(uwlinca_debug) {fprintf(stderr,"LOCATION==%d(fast-y, not bottom-up)\n", location); }
        }
    } else {  // fast-X, cca data
        if ( strcmp(uwlinca_configuration->seek_axis, "fast-x") == 0 ||
                     strcmp(uwlinca_configuration->seek_axis, "fast-X") == 0 ) { // fast-x,  uwlinca 
            if(strcmp(uwlinca_configuration->seek_direction, "bottom-up") == 0) { 
                    location = ((long)z * uwlinca_configuration->nx * uwlinca_configuration->ny) + (y * uwlinca_configuration->nx) + x;
if(uwlinca_debug) {fprintf(stderr,"LOCATION==%d(fast-x, bottom-up)\n", location); }
                } else { // bottom-up
                        location = ((long)(uwlinca_configuration->nz - z) * uwlinca_configuration->nx * uwlinca_configuration->ny) + (y * uwlinca_configuration->nx) + x;
if(uwlinca_debug) {fprintf(stderr,"LOCATION==%d(fast-x, not bottom-up)\n", location); }
            }
        }
    }

    // Check our loaded components of the model.
    if (uwlinca_velocity_model->vp_status == 2) {
        // Read from memory.
        ptr = (float *)uwlinca_velocity_model->vp;
        data->vp = ptr[location];
    } else if (uwlinca_velocity_model->vp_status == 1) {
        // Read from file.
        fp = (FILE *)uwlinca_velocity_model->vp;
        fseek(fp, location * sizeof(float), SEEK_SET);
        float temp;
        fread(&(temp), sizeof(float), 1, fp);
if(uwlinca_debug) {fprintf(stderr,"     FOUND : vp %f\n", temp); }
        data->vp=temp;
    }

    /* Calculate vs */
    if (data[p].vp > 0.0) { data[p].vs=uwlinca_calculate_vs(data[p].vp); }
    /* Calculate density */
    if (data[p].vp > 0.0) { data[[].rho=uwlinca_calculate_density(data[p].vp); }
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void uwlinca_trilinear_interpolation(double x_percent, double y_percent, double z_percent,
                             uwlinca_properties_t *eight_points, uwlinca_properties_t *ret_properties) {
    uwlinca_properties_t *temp_array = calloc(2, sizeof(uwlinca_properties_t));
    uwlinca_properties_t *four_points = eight_points;

    uwlinca_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

    // Now advance the pointer four "cvms5_properties_t" spaces.
    four_points += 4;

    // Another interpolation.
    uwlinca_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

    // Now linearly interpolate between the two.
    uwlinca_linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

    free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void uwlinca_bilinear_interpolation(double x_percent, double y_percent, uwlinca_properties_t *four_points, uwlinca_properties_t *ret_properties) {
    uwlinca_properties_t *temp_array = calloc(2, sizeof(uwlinca_properties_t));
    uwlinca_linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
    uwlinca_linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
    uwlinca_linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
    free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void uwlinca_linear_interpolation(double percent, uwlinca_properties_t *x0, uwlinca_properties_t *x1, uwlinca_properties_t *ret_properties) {
    ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
    ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
    ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
    ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
    ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int uwlinca_finalize() {

    proj_destroy(uwlinca_geo2utm);
    uwlinca_geo2utm = NULL;

    if (uwlinca_velocity_model) free(uwlinca_velocity_model);
    if (uwlinca_configuration) free(uwlinca_configuration);

    return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int uwlinca_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(uwlinca_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, uwlinca_version_string, verlen);
  return 0;
}

/**
 * Returns the model config information.
 *
 * @param key Config key string to return.
 * @param sz number of config terms.
 * @return Zero
 */
int uwlinca_config(char **config, int *sz)
{
  int len=strlen(uwlinca_config_string);
  if(len > 0) {
    *config=uwlinca_config_string;
    *sz=uwlinca_config_sz;
    return SUCCESS;
  }
  return FAIL;
}


/**
 * Reads the uwlinca_configuration file describing the various properties of CVM-S5 and populates
 * the uwlinca_configuration struct. This assumes uwlinca_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The uwlinca_configuration file location on disk to read.
 * @param config The uwlinca_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int uwlinca_read_configuration(char *file, uwlinca_configuration_t *config) {
    FILE *fp = fopen(file, "r");
    char key[40];
    char value[80];
    char line_holder[128];

    // If our file pointer is null, an error has occurred. Return fail.
    if (fp == NULL) {
        uwlinca_print_error("Could not open the uwlinca_configuration file.");
        return FAIL;
    }

    // Read the lines in the uwlinca_configuration file.
    while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
        if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
            sscanf(line_holder, "%s = %s", key, value);

            // Which variable are we editing?
            if (strcmp(key, "utm_zone") == 0) config->utm_zone = atoi(value);
            if (strcmp(key, "model_dir") == 0) sprintf(config->model_dir, "%s", value);
            if (strcmp(key, "nx") == 0) config->nx = atoi(value);
            if (strcmp(key, "ny") == 0) config->ny = atoi(value);
            if (strcmp(key, "nz") == 0) config->nz = atoi(value);
            if (strcmp(key, "depth") == 0) config->depth = atof(value);
            if (strcmp(key, "top_left_corner_e") == 0) config->top_left_corner_e = atof(value);
            if (strcmp(key, "top_left_corner_n") == 0) config->top_left_corner_n = atof(value);
            if (strcmp(key, "top_right_corner_e") == 0) config->top_right_corner_e = atof(value);
            if (strcmp(key, "top_right_corner_n") == 0) config->top_right_corner_n = atof(value);
            if (strcmp(key, "bottom_left_corner_e") == 0) config->bottom_left_corner_e = atof(value);
            if (strcmp(key, "bottom_left_corner_n") == 0) config->bottom_left_corner_n = atof(value);
            if (strcmp(key, "bottom_right_corner_e") == 0) config->bottom_right_corner_e = atof(value);
            if (strcmp(key, "bottom_right_corner_n") == 0) config->bottom_right_corner_n = atof(value);
            if (strcmp(key, "depth_interval") == 0) config->depth_interval = atof(value);
            if (strcmp(key, "seek_axis") == 0) sprintf(config->seek_axis, "%s", value);
            if (strcmp(key, "seek_direction") == 0) sprintf(config->seek_direction, "%s", value);
            if (strcmp(key, "interpolation") == 0) { 
                config->interpolation=0;
                if (strcmp(value,"on") == 0) config->interpolation=1;
            }
        }
    }

    // Have we set up all uwlinca_configuration parameters?
    if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' || 
        config->seek_direction[0] == '\0' || config->seek_axis[0] == '\0' ||
        config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
        config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
        config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
        config->depth_interval == 0) {
        uwlinca_print_error("One uwlinca_configuration parameter not specified. Please check your uwlinca_configuration file.");
        return FAIL;
    }

    fclose(fp);

    return SUCCESS;
}

/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void uwlinca_print_error(char *err) {
    fprintf(stderr, "An error has occurred while executing uwlinca. The error was: %s\n",err);
    fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
    fprintf(stderr, "about the computer you are running uwlinca on (Linux, Mac, etc.).\n");
}

/**
 * Check if the data is too big to be loaded internally (exceed maximum
 * allowable by a INT variable)
 *
 */
static int too_big() {
        long max_size= (long) (uwlinca_configuration->nx) * uwlinca_configuration->ny * uwlinca_configuration->nz;
        long delta= max_size - INT_MAX;

    if( delta > 0) {
        return 1;
        } else {
        return 0;
        }
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int uwlinca_try_reading_model(uwlinca_model_t *model) {
    double base_malloc = uwlinca_configuration->nx * uwlinca_configuration->ny * uwlinca_configuration->nz * sizeof(float);
    int file_count = 0;
    int all_read_to_memory =0;
    char current_file[128];
    FILE *fp;

    // Let's see what data we actually have.
    sprintf(current_file, "%s/vp.dat", uwlinca_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
            model->vp = malloc(base_malloc);
            if (model->vp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vp, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->vp_status = 2;
            } else {
              model->vp = fopen(current_file, "rb");
              model->vp_status = 1;
            }
        } else {
            model->vp = fopen(current_file, "rb");
            model->vp_status = 1;
                }
        file_count++;
    }

    if (file_count == 0)
        return FAIL;
    else if (file_count > 0 && all_read_to_memory != file_count)
        return SUCCESS;
    else
        return 2;
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls uwlinca_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
    return uwlinca_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls uwlinca_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(uwlinca_point_t *points, uwlinca_properties_t *data, int numpoints) {
    return uwlinca_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls uwlinca_finalize.
 *
 * @return Success
 */
int model_finalize() {
    return uwlinca_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls uwlinca_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
    return uwlinca_version(ver, len);
}

/**
 * Version function loaded and called by the UCVM library. Calls uwlinca_config.
 *
 * @param config Config string to return.
 * @param sz number of config terms
 * @return Zero
 */
int model_config(char **config, int *sz) {
    return uwlinca_config(config, sz);
}


int (*get_model_init())(const char *, const char *) {
        return &uwlinca_init;
}
int (*get_model_query())(uwlinca_point_t *, uwlinca_properties_t *, int) {
         return &uwlinca_query;
}
int (*get_model_finalize())() {
         return &uwlinca_finalize;
}
int (*get_model_version())(char *, int) {
         return &uwlinca_version;
}
int (*get_model_config())(char **, int*) {
    return &uwlinca_config;
}



#endif
