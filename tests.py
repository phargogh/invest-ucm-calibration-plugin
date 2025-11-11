import shutil
import tempfile
import unittest

from invest_ucm_calibration import plugin


class UCMCalibrationPluginTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workspace)

    def test_execute(self):
        plugin.execute({})
