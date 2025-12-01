import logging
import os
import pprint

from natcap.invest import gettext
from natcap.invest import spec
from natcap.invest import utils
from natcap.invest import validation
from natcap.invest.unit_registry import u

LOGGER = logging.getLogger(__name__)

MODEL_SPEC = spec.ModelSpec(
    model_id="invest-ucm-calibration",
    model_title="Urban Cooling Model Calibration",
    userguide=(
        "https://invest-ucm-calibration.readthedocs.io/en/latest/usage.html"),
    module_name=__name__,
    input_field_order=[
        ['workspace_dir'],
        ['lulc_raster_path', 'aoi_vector_path'],
        ['cc_method', 'ref_eto_table'],
        ['t_rasters_table', 't_stations', 'uhi_max'],
        ['metric', 'stepsize', 'exclude_zero_kernel_dist'],
        ['num_steps', 'num_update_logs'],
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
        spec.OptionStringInput(
            id="cc_method",
            name=gettext("cooling capacity calculation method"),
            about=gettext("The air temperature predictor method to use."),
            options=[
                spec.Option(
                    key="factors",
                    about=(
                        "Use the weighted shade, albedo, and ETI factors as a temperature"
                        " predictor (for daytime temperatures).")),
                spec.Option(
                    key="intensity",
                    about=(
                        "Use building intensity as a temperature predictor (for nighttime"
                        " temperatures)."))
            ]
        ),
        spec.CSVInput(
            id="ref_eto_table",
            name=gettext("Reference evapotranspiration rasters"),
            about=gettext("Table of reference evapotranspiration rasters"),
            index_col="eto_path",
            columns=[
                spec.SingleBandRasterInput(
                    id="eto_path",
                    name=gettext('Reference evapotranspiration raster'),
                    about=gettext(
                        "The path to a raster of reference "
                        "evapotranspiration values."),
                    data_type=float,
                    units=None,
                ),
            ]
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
        spec.NumberInput(
            id="uhi_max",
            name=gettext("UHI effect"),
            required=False,  # can be interred from temp rasters
            about=gettext(
                "The magnitude of the urban heat island effect, i.e., the difference"
                " between the rural reference temperature and the maximum temperature"
                " observed in the city. This model is designed for cases where"
                " UHI is positive, meaning the urban air temperature is greater"
                " than the rural reference temperature."
            ),
            units=u.degree_Celsius,
            expression="value >= 0",
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
                "Defaults to TODO"),
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
        spec.AOI.model_copy(update=dict(id="aoi_vector_path")),
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
    calibrator_args = {}
    calibrator_args.update({
        'dst_filepath': os.path.join(
            args['workspace_dir'], 'calibration-results.json'),
    })
    pprint.pprint(args)
    pprint.pprint(calibrator_args)

    #if not args['uhi_max']:
    #    # Calculate from the max/min observed temps, for each station/date
    #    pass


@validation.invest_validator
def validate(args, limit_to=None):
    return validation.validate(args, MODEL_SPEC)
