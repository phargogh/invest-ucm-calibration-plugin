from natcap.invest import gettext
from natcap.invest import spec
from natcap.invest import utils
from natcap.invest import validation
from natcap.invest.unit_registry import u

MODEL_SPEC = spec.ModelSpec(
    model_id="invest-ucm-calibration",
    model_title="Urban Cooling Model Calibration",
    userguide=(
        "https://invest-ucm-calibration.readthedocs.io/en/latest/usage.html"),
    input_field_order=[
    ],
    inputs=[
        spec.DirectoryInput(
            id="workspace_dir",
            name=gettext("workspace"),
            about=gettext(
                "The folder where all the model's output files will be "
                "written. If this folder does not exist, it will be created. "
                "If data already exists in the folder, it will be "
                "overwritten."),
            contents={},
            must_exist=False,
            permissions="rwx"
        ),
    ],
    outputs=[
    ],
)


def execute(args):
    pass


@validation.invest_validator
def validate(args):
    return validation.validate(args, MODEL_SPEC)
