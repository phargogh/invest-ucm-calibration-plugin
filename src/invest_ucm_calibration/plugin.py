import logging
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
    pprint.pprint(args)


@validation.invest_validator
def validate(args, limit_to=None):
    return validation.validate(args, MODEL_SPEC)
