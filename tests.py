import os
import shutil
import subprocess
import tempfile
import unittest

from invest_ucm_calibration import plugin

CWD = os.path.dirname(__file__)
DATA = os.path.join(CWD, 'data')
DATA_ZIP = 'https://github.com/martibosch/invest-ucm-calibration/archive/refs/tags/v0.6.0.zip'

if not os.path.exists(DATA):
    print(f"Downloading data directory from {DATA_ZIP}")
    subprocess.run(['wget', DATA_ZIP, -O, 'data.zip'])
    subprocess.run(['unzip', 'tests/data', '-d', DATA])
    os.remove('data.zip')


class UCMCalibrationPluginTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workspace)

    def test_execute(self):
        args = {
            'workspace_dir': self.workspace_dir,
            'lulc_raster_path': 'foo',
            'cc_method': 'intensity',
            'aoi_vector_filepath': '',
        }
        plugin.execute(args)

    def test_validate(self):
        args = {
            'workspace_dir': self.workspace_dir,
            'lulc_raster_path': '',
            'cc_method': 'intensity',
            'aoi_vector_path': ''
        }
        self.assertEqual([], plugin.validate(args))
