import argparse
import os
import re
import textwrap

import numpy
import pygeoprocessing
from natcap.invest import spec
from osgeo import gdal
from osgeo import ogr


def render_raster(raster_path, workspace_path):
    stem, ext = os.path.splitext(os.path.basename(raster_path))

    raster_info = pygeoprocessing.get_raster_info(raster_path)
    if raster_info['n_bands'] != 1:
        raise NotImplementedError("Can't yet handle multiband rasters")

    numpy_dtype = raster_info['numpy_type']

    if numpy.issubdtype(numpy_dtype, numpy.floating):
        data_type = 'float'
    else:
        data_type = 'int'

    return textwrap.dedent(
        f"""\
        spec.SingleBandRasterOutput(
            id="{stem}",
            path="{workspace_path}",
            about=gettext(""),  # FIXME
            data_type={data_type},
            units=None,  # FIXME
        ),
        """)


def render_vector(vector_path, workspace_path):
    stem, ext = os.path.splitext(os.path.basename(vector_path))
    try:
        vector = gdal.OpenEx(vector_path)
        layer = vector.GetLayer()

        fields = []
        for field in layer.schema:
            if field.GetType() in {ogr.OFTInteger, ogr.OFTInteger64}:
                spec_classname = "IntegerOutput"
            elif field.GetType() in {ogr.OFTReal}:
                spec_classname = "NumberOutput"
            else:  # assume it's a string
                spec_classname = "StringOutput"
            fields.append(textwrap.dedent(
                f"""\
                spec.{spec_classname}(  # TODO: Confirm correct type
                    id="{field.GetName()}",
                    about=gettext(""),  # FIXME
                    units=None,  # FIXME
                ),
                """))

        # Compute the geometry name
        geom_types = set()
        ogr_geom_types = {}
        for attrname in filter(lambda a: a.startswith('wkb'), dir(ogr)):
            new_attrname = re.sub('^wkb', '', attrname)
            new_attrname = re.sub('[A-Z0-9]+$', '', new_attrname)
            ogr_geom_types[getattr(ogr, attrname)] = new_attrname.upper()

        for geom_type, feature_count in layer.GetGeometryTypes().items():
            geom_types.add(ogr_geom_types[geom_type])
    finally:
        layer = None
        vector = None

    value = textwrap.dedent(
        f"""\
        spec.VectorOutput(
            id="{stem}",
            path="{workspace_path}",
            about=gettext(""),  # FIXME
            geometry_types={sorted(geom_types)},
            fields=[
            %s
            ],
        ),
        """)

    # Textwrap.indent inside of a dedent call was giving me weird beahvior, but
    # this works well.
    value = value % textwrap.indent(''.join(fields), prefix='    ')
    return value


# build a set of known output spec classes
def main():
    parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__),
        description=(
            'Build an output spec given a workspace produced by an '
            'InVEST-like model, such as a plugin.')
    )
    parser.add_argument('workspace', help='The workspace to analyze')

    # TODO: optionally write to an output file.
    args = parser.parse_args()

    dirname = args.workspace
    # ASSUMPTION: `dirname` is the workspace directory itself; everything will
    # be relative to this workspace.
    dirname = os.path.expanduser(dirname)  # in case ~ is used

    extensions_and_templates = {
        '.tif': render_raster,
        '.shp': render_vector,
        '.db': lambda x, y: "spec.TASKGRAPH_CACHE,\n",
    }
    extensions_to_exclude = {'.prj', '.dbf', '.shx'}

    output_spec = []
    for root, dirs, files in os.walk(dirname):
        # root is the current directory being visited
        # dirs is the list of subdirectories in root
        # files is the list of files in root

        # iterate over files first, then iterate over directories later.
        for filename in files:
            filepath = os.path.join(root, filename)
            relpath = os.path.relpath(filepath, dirname)
            stem, ext = os.path.splitext(filepath)
            if ext in extensions_to_exclude:
                print(f"Excluding {ext} file {relpath}")
            elif ext not in extensions_and_templates:
                print(f"Don't know how to handle {ext} for file {relpath}")
            else:
                templated_string = extensions_and_templates[ext](
                    filepath, relpath)
                output_spec.append(templated_string)

    print(''.join(output_spec))


if __name__ == '__main__':
    main()
