import glob
import os
import shutil
import subprocess
import tempfile
import unittest

import pygeoprocessing
from invest_ucm_calibration import plugin
from osgeo import gdal
from osgeo import osr

CWD = os.path.dirname(__file__)
DATA = os.path.join(CWD, 'data')
DATA_ZIP = 'https://github.com/martibosch/invest-ucm-calibration/archive/refs/tags/v0.6.0.zip'

if not os.path.exists(DATA):
    print(f"Downloading data directory from {DATA_ZIP}")
    LOCAL_ZIP = os.path.join(CWD, 'data.zip')
    subprocess.run(['wget', DATA_ZIP, '-O', LOCAL_ZIP], check=True)
    subprocess.run(['unzip', '-j', LOCAL_ZIP,
                    'invest-ucm-calibration-0.6.0/tests/data/*', '-d', DATA],
                   check=True)
    assert len(os.listdir(DATA)), 'Data did not expand correctly'
    os.remove(LOCAL_ZIP)

    # LULC is shipped with a geographic coordinate system, which isn't allowed
    # by the model.  Warp to a local projection.
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(32632)  # WGS84 / UTM zone 32 N
    ten_km = 1000*10
    pygeoprocessing.warp_raster(
        os.path.join(DATA, 'lulc.tif'), [ten_km, -ten_km],
        os.path.join(DATA, 'lulc_utm32n.tif'), resample_method='near',
        target_projection_wkt=target_srs.ExportToWkt())


# Parameters here match what we use in UCM, not the ones that the calibration
# tool uses.  The plugin will handle the mapping.
TEST_KWARGS = {
    # standard args for UCM
    'lulc_raster_path': os.path.join(DATA, 'lulc_utm32n.tif'),
    'biophysical_table_path': os.path.join(DATA, 'biophysical-table.csv'),
    'aoi_vector_path': os.path.join(DATA, 'aoi.gpkg'),
    'cc_method': 'factors',  # or 'intensity'

    # Calibration tool takes modified type (list)
    'ref_eto_raster_path': glob.glob(os.path.join(DATA, 'ref_et*.tif'))[0],
    't_ref': glob.glob(os.path.join(DATA, 'T*.tif'))[0],  # CALTOOL: takes list

    # nonstandard args, for the calibration tool
    'station_t_filepath': os.path.join(DATA, 'station-t.csv'),
    'station_locations': os.path.join(DATA, 'station-locations.csv'),
    'station_t_one_day_filepath': os.path.join(DATA, 'station-t-one-day.csv'),
    'num_steps': 2,
    'num_update_logs': 2,

}


class UCMCalibrationPluginTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workspace)

    def test_execute(self):
        args = {
            'workspace_dir': self.workspace,
            'lulc_raster_path': 'foo',
            'cc_method': 'intensity',
            'aoi_vector_filepath': '',
        }
        plugin.execute(args)

    def test_validate(self):
        args = TEST_KWARGS.copy()
        args['workspace_dir'] = self.workspace
        self.assertEqual([], plugin.validate(args))
