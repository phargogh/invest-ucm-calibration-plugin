import collections
import logging
import os
import pprint
import re

from invest_ucm_calibration import settings as ucm_cal_defaults
from invest_ucm_calibration.cli import main as ucm_cal_main
from natcap.invest import gettext
from natcap.invest import spec
from natcap.invest import validation
from natcap.invest.unit_registry import u
from osgeo import gdal

# TODO: should I be using the calibration tool CLI or the object?
#       CLI controls logging, warnings

LOGGER = logging.getLogger(__name__)

MODEL_SPEC = spec.ModelSpec(
    model_id="invest-ucm-calibration",
    model_title="Urban Cooling Model Calibration",
    userguide=(
        "https://invest-ucm-calibration.readthedocs.io/en/latest/usage.html"),
    module_name=__name__,
    input_field_order=[
        ['workspace_dir'],
        ['lulc_raster_path', 'biophysical_table_path', 'aoi_vector_path'],
        ['cc_method', 'ref_eto_table'],
        ['t_refs', 't_rasters_table', 't_stations', 'uhi_maxs'],
        ['metric', 'stepsize', 'exclude_zero_kernel_dist'],
        ['num_steps', 'num_update_logs', 'initial_solution'],
    ],
    inputs=[
        spec.WORKSPACE,
        spec.SingleBandRasterInput(
            id="lulc_raster_path",
            name=gettext("Land use/land cover"),
            about=gettext(
                "Map of land use / land cover codes. Each land use/land "
                "cover type must be assigned a unique integer code."
            ),
            required=True,
            data_type=int,  # TODO is there a preferred type for this?
            units=None,
            projected=True,
            projection_units=u.meter,
        ),
        spec.CSVInput(
            id="biophysical_table_path",
            name=gettext("biophysical table"),
            about=gettext(
                "A table mapping each LULC code to biophysical data for that "
                "LULC class. All values in the LULC raster must have "
                "corresponding entries in this table."
            ),
            columns=[
                spec.LULC_TABLE_COLUMN,
                spec.NumberInput(
                    id="kc",
                    about=gettext("Crop coefficient for this LULC class."),
                    units=u.none
                ),
                spec.BooleanInput(
                    id="green_area",
                    about=gettext(
                        "Enter 1 to indicate that the LULC is considered a "
                        "green area. Enter 0 to indicate that the LULC is "
                        "not considered a green area."
                    )
                ),
                spec.RatioInput(
                    id="shade",
                    about=(
                        "The proportion of area in this LULC class that is "
                        "covered by tree canopy at least 2 meters high. "
                        "Required if the 'factors' option is selected for the "
                        "Cooling Capacity Calculation Method."
                    ),
                    required="cc_method == 'factors'",
                    units=None
                ),
                spec.RatioInput(
                    id="albedo",
                    about=(
                        "The proportion of solar radiation that is directly "
                        "reflected by this LULC class. Required if the "
                        "'factors' option is selected for the Cooling "
                        "Capacity Calculation Method."
                    ),
                    required="cc_method == 'factors'",
                    units=None
                ),
                spec.RatioInput(
                    id="building_intensity",
                    about=(
                        "The ratio of building floor area to footprint area, "
                        "with all values in this column normalized between 0 "
                        "and 1. Required if the 'intensity' option is "
                        "selected for the Cooling Capacity Calculation Method."
                    ),
                    required="cc_method == 'intensity'",
                    units=None
                )
            ],
            index_col="lucode"
        ),
        spec.OptionStringInput(
            id="cc_method",
            name=gettext("cooling capacity calculation method"),
            about=gettext("The air temperature predictor method to use."),
            options=[
                spec.Option(
                    key="factors",
                    about=(
                        "Use the weighted shade, albedo, and ETI factors as a "
                        "temperature predictor (for daytime temperatures).")),
                spec.Option(
                    key="intensity",
                    about=(
                        "Use building intensity as a temperature predictor "
                        "(for nighttime temperatures)."))
            ]
        ),
        spec.CSVInput(
            id="ref_eto_table",
            name=gettext("Reference evapotranspiration rasters"),
            about=gettext("Table of reference evapotranspiration rasters"),
            columns=[
                spec.SingleBandRasterInput(
                    id="eto_path",
                    name=gettext('Reference evapotranspiration raster'),
                    about=gettext(
                        "The path to a raster of reference "
                        "evapotranspiration values."),
                    data_type=str,
                    units=None,
                ),
            ]
        ),
        # If an AOI is not provided, the calibration tool uses the bounds of
        # the LULC instead.
        spec.AOI.model_copy(update=dict(id="aoi_vector_path", optional=True)),

        spec.StringInput(
            id="t_refs",
            name=gettext("Reference Air Temperatures"),
            about=gettext(
                "Reference air temperature or temperatures.  If multiple "
                "temperatures are provided, they must be comma-separated."
            ),
            regexp="[0-9., ]+",
            expression="float(v.strip()) for v in value.split(',')",
            # If not provided, the calibration tool uses the minimum observed
            # temperature, either from raster or station measurements,
            # depending on which sets of inputs are provided.
            optional=True,
            units=u.degree_Celsius,
        ),
        spec.NumberInput(
            id="uhi_maxs",
            name=gettext("UHI effects"),
            required=False,  # can be interred from temp rasters
            about=gettext(
                "The magnitude of the urban heat island effect, i.e., the "
                "difference between the rural reference temperature and the "
                "maximum temperature observed in the city. This model is "
                "designed for cases where UHI is positive, meaning the urban "
                "air temperature is greater than the rural reference "
                "temperature.  If not provided, the difference between the "
                "minimum and maximum observed temperature will be used, "
                "calculated from input temperature rasters or station "
                "measurements, for each respective date if calibrating for "
                "multiple dates.  If providing multiple UHI values, such as "
                "for multiple dates, UHI values must be comma-separated."
            ),
            units=u.degree_Celsius,
            regexp="[0-9., ]+",
            expression="all(float(v.strip()) > 0 for v in value.split(','))",
        ),
        spec.StringInput(
            id="initial_solution",
            name=gettext("Initial Solution"),
            about=gettext(
                "Sequence with the parameter values used as initial "
                "solution, which can either be of the form "
                "(t_air_average_radius, green_area_cooling_distance, "
                "cc_weight_shade, cc_weight_albedo, cc_weight_eti) when "
                "`cc_method` is 'factors', or (t_air_average_radius, "
                "green_area_cooling_distance) when `cc_method` is "
                "'intensity'. If not provided, the default values of the "
                "urban cooling model will be used."),
            regexp="[0-9., ]+",
            expression="float(v.strip()) for v in value.split(',')",
            required=False,
        ),
        spec.CSVInput(
            id="t_rasters_table",
            name=gettext("Table of temperature rasters"),
            about=gettext(
                "Table of temperature rasters and observation dates."),
            index_col="t_raster_date",
            required="not t_stations",
            columns=[
                spec.StringInput(
                    id="t_raster_date",
                    name=gettext("Date of temperature observation."),
                    about=gettext(
                        "The date of temperature observation, in the form "
                        "DD-MM-YYYY"),
                    regexp="[0-3][0-9]-[0-1][0-9]-[1-2][0-9][0-9][0-9]",
                    expression="str(value)",
                ),
                spec.SingleBandRasterInput(
                    id="t_raster_path",
                    name=gettext("Temperature raster path"),
                    about=gettext(
                        "The path to a raster of air temperatures observed "
                        "on the associated date."),
                    data_type=float,
                    units=u.degree_Celsius,
                ),
            ]
        ),
        spec.VectorInput(
            id="t_stations",
            name=gettext("Temperature stations vector"),
            about=gettext(
                "Spatial vector of temperature observations from stations."),
            index_col="",
            required="not t_rasters_table",
            geometry_types={'POINT'},
            fields=[
                spec.StringInput(
                    id="name",
                    required=False,
                    name=gettext("Station name"),
                    about=gettext("The name of the monitoring station."),
                ),
                spec.NumberInput(
                    id="[DATE]",
                    required=True,
                    name=gettext("Temperature recording"),
                    about=gettext("The date the temperature was recorded."),
                    units=u.degree_Celsius,
                    exression="float(value)",
                ),
            ],
            projected=True,  # will this be necessary?
        ),
        spec.OptionStringInput(
            id="metric",
            name=gettext("Metric to optimize"),
            about=gettext("Target metric to optimize in the calibration."),
            options=[
                spec.Option(
                    key="R2",
                    about=(
                        "Maximize the R squared metric during calibration.")),
                spec.Option(
                    key="MAE",
                    about=(
                        "Minimize the Mean Absolute Error during "
                        "calibration.")),
                spec.Option(
                    key="RMSE",
                    about=(
                        "Minimize the root mean squared error during "
                        "calibration.")),
            ],
        ),
        spec.NumberInput(
            id="stepsize",
            name="Calibration step size",
            about=(
                # Copied from the calibration tool.
                "Step size in terms of the fraction of each parameter when "
                "looking to select a neighbor solution for the following "
                "iteration. The neighbor will be randomly drawn from an "
                "uniform distribution in the [param - stepsize * param, "
                "param + stepsize * param] range. For example, with a step "
                "size of 0.3 and a `t_air_average_radius` of 500 at a given"
                "iteration, the solution for the next iteration will be "
                "uniformly sampled from the [350, 650] range."),
            units=None,
            expression="float(value) > 0",
        ),
        spec.BooleanInput(
            id="exclude_zero_kernel_dist",
            name=gettext("Exclude zero-size kernels"),
            about=gettext(
                "Whether the calibration should consider parameters that "
                "lead to decay functions with a kernel distance of zero "
                "pixels (i.e., `t_air_average_radius` or "
                "`green_area_cooling_distance` lower than half the LULC pixel "
                "resolution)."),
        ),
        spec.IntegerInput(
            id="num_steps",
            name=gettext("Number of Calibration Steps"),
            about=gettext(
                "Number of steps in the simulated annealing procedure. "
                f"Defaults to {ucm_cal_defaults.DEFAULT_NUM_STEPS}"),
            expression="int(value) > 0",
        ),
        spec.IntegerInput(
            id="num_update_logs",
            name=gettext("Number of updates logged"),
            about=gettext(
                "The number of updates to log during simulated annealing. "
                "If this number is the same as the number of steps, then "
                "each iteration will be logged"),
            expression="int(value) > 0",
        ),

        # TODO ref_et_raster_filepaths
        # T_REFs - numeric or iterable of numbers
    ],
    outputs=[
    ],
)

# Checkboxes
#   * Calibrate with dates?  Yes/No
#   * Calibrate against temperature maps?     | One or the other
#   * Calibrate against station measurements? | One or the other

# New table: dates
# Columns:
#   * date (probably use an ISO-compatible datetime) (required if calib w/dates)
#       --> if we aren't calibrating with dates, then we just pass None to the
#           calibration tool and let the tool do its thing
#   * t_ref (req)
#   * uhi_max (optional)
#   * t_raster_filepaths (required if calibrating against temp maps)
#   * ref_et_raster_filepaths (implied by source code to be associated with
#     dates, but perhaps dates might be arange, as below?)

# Notes:
#   * user-provided dates are not used if a station filepath is provided
#       * dates are taken from the station table instead
#   * if the user does not provide dates, we assume some integer dates
#       (numpy.arange(len(ref_et_raster_filepaths)))


def execute(args):
    calibrator_defaults = {}
    for attrname in dir(ucm_cal_defaults):
        if not attrname.startswith('DEFAULT_'):
            continue
        value = getattr(ucm_cal_defaults, attrname)
        calibrator_defaults[attrname.lower().replace('default_', '')] = value

    calibrator_args = {}
    calibrator_args.update({
        'dst_filepath': os.path.join(
            args['workspace_dir'], 'calibration-results.json'),
    })
    pprint.pprint(args)
    pprint.pprint(calibrator_args)
    pprint.pprint(calibrator_defaults)

    # Some arguments can be copied directly over to the calibrator args
    for plugin_key, calibrator_key in [
            ('lulc_raster_path', 'lulc_raster_filepath'),
            ('biophysical_table_path', 'biophysical_table_filepath'),
            ('cc_method', 'cc_method'),]:
        calibrator_args[calibrator_key] = args[plugin_key]
    # Translate the Ref ET0 vector to the lists that the calibrator expects.
    # TODO: use InVEST's table loading
    et0_df = MODEL_SPEC.get_input('ref_eto_table').get_validated_dataframe(
        args['ref_eto_table'])
    # The calibrator CLI expects a comma-separated string list of paths
    calibrator_args['ref_et_raster_filepaths'] = ','.join(list(
        et0_df['eto_path'].values))

    # Translate t_stations vector to the CSV expected by the plugin
    # TODO: use file registry
    # TODO: use taskgraph
    stations_loc_csv = os.path.join(args['workspace_dir'], 'station-loc.csv')
    stations_temps_csv = os.path.join(
        args['workspace_dir'], 'station-temps.csv')
    _t_stations_vector_to_csv(
        args['t_stations'], stations_loc_csv, stations_temps_csv)
    calibrator_args['station_locations_filepath'] = stations_loc_csv
    calibrator_args['station_t_filepath'] = stations_temps_csv

    # If not provided, default model args will be used.
    # TODO: How does the calibration tool actually handle initial_solution??
    #       Do we need to provide an initial solution or will the calibrator do
    #       it for us?
    for key in ('initial_solution', 't_refs', 'uhi_maxs'):
        if key in args and args[key]:
            value = args[key]
            if isinstance(value, str):
                value = value.split(',')
            elif isinstance(value, (int, float)):
                value = [value]
            # float() will chomp leading/trailing whitespace if present
            calibrator_args[key] = [float(v) for v in value]

    ucm_cal_main.cli(**calibrator_args)

    # if not args['uhi_max']:
    #     # Calculate from the max/min observed temps, for each station/date
    #     pass


def _t_stations_vector_to_csv(
        t_stations_vector_path, t_stations_locations_csv_path,
        t_stations_temps_csv_path):
    vector = gdal.OpenEx(t_stations_vector_path)
    layer = vector.GetLayer()

    # TODO: do I need to transform the points from local projection to WGS84?

    stations = {}  # name: (x, y)
    temps = collections.defaultdict(dict)  # date: name: temp
    for feature in layer:
        geom = feature.GetGeometryRef()
        name = feature.GetField('name')
        stations[name] = (geom.GetX(), geom.GetY())

        for datefieldname, field_value in feature.items().items():
            # Match various forms of date, e.g. YYYY-MM-DD, DD-MM-YYYY
            if not re.match('[0-9-]+', datefieldname):
                continue
            temps[datefieldname][name] = field_value

    with open(t_stations_locations_csv_path, 'w') as station_locations:
        station_locations.write(",x,y\n")
        for st_name, (x, y) in stations.items():
            station_locations.write(f"{st_name},{x},{y}\n")

    with open(t_stations_temps_csv_path, 'w') as station_temps:
        station_names = list(stations.keys())
        names_header = ','.join(station_names)
        station_temps.write(f',{names_header}\n')
        for date, station_data in temps.items():
            temps_data = ','.join(str(station_data[name]) for name in station_names)
            station_temps.write(f'{date},{temps_data}\n')


@validation.invest_validator
def validate(args, limit_to=None):
    return validation.validate(args, MODEL_SPEC)
