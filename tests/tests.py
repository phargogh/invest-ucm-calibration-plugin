import glob
import os
import random
import shutil
import subprocess
import tempfile
import textwrap
import unittest
import urllib.request

import numpy
import pandas
import pandas.testing
import pygeoprocessing
import scipy.ndimage
import shapely.geometry
from hypothesis import given
from hypothesis import strategies as st
from invest_ucm_calibration_plugin import plugin
from natcap.invest import datastack
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from shapely.geometry import Point
from unittest_parametrize import parametrize
from unittest_parametrize import ParametrizedTestCase

CWD = os.path.dirname(__file__)
DATA = os.path.join(CWD, 'data')
DATA_ZIP = (
    'https://github.com/martibosch/invest-ucm-calibration/archive/'
    'refs/tags/v0.6.0.zip')
INVEST_DATA = os.path.join(CWD, 'invest-data')
INVEST_DATA_ZIP = (
    'https://storage.googleapis.com/releases.naturalcapitalproject.org/'
    'invest/3.17.2/data/UrbanCoolingModel.zip')


def _randomize_pts_in_bbox(bbox, n_pts, seed=None):
    if seed:
        random.seed(seed)
    minx, miny, maxx, maxy = bbox

    return [shapely.geometry.Point([
        random.uniform(minx, maxx), random.uniform(miny, maxy)]) for i in
        range(n_pts)]


class UCMCalibrationPluginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(DATA):
            LOCAL_ZIP = os.path.join(CWD, 'data.zip')
            subprocess.run(['wget', DATA_ZIP, '-O', LOCAL_ZIP], check=True)
            subprocess.run(['unzip', '-j', LOCAL_ZIP,
                            'invest-ucm-calibration-0.6.0/tests/data/*',
                            '-d', DATA],
                           check=True)
            assert len(os.listdir(DATA)), 'Data did not expand correctly'
            os.remove(LOCAL_ZIP)

            # LULC is shipped with a geographic coordinate system, which isn't
            # allowed by the model.  Warp to a local projection.
            target_srs = osr.SpatialReference()
            target_srs.ImportFromEPSG(32632)  # WGS84 / UTM zone 32 N
            ten_km = 1000*10
            pygeoprocessing.warp_raster(
                os.path.join(DATA, 'lulc.tif'), [ten_km, -ten_km],
                os.path.join(DATA, 'lulc_utm32n.tif'), resample_method='near',
                target_projection_wkt=target_srs.ExportToWkt())

        if not os.path.exists(INVEST_DATA):
            with tempfile.TemporaryDirectory() as temp_dir:
                zipfile_url = INVEST_DATA_ZIP
                target_zipfile = os.path.join(temp_dir, 'zipfile.zip')
                print(f"Downloading data directory from {zipfile_url}")
                urllib.request.urlretrieve(zipfile_url, target_zipfile)
                shutil.unpack_archive(target_zipfile, INVEST_DATA)

    def setUp(self):
        self.workspace = tempfile.mkdtemp()

        # Parameters here match what we use in UCM, not the ones that the
        # calibration tool uses.  The plugin will handle the mapping.
        self.base_kwargs = {
            # standard args for UCM
            'lulc_raster_path': os.path.join(DATA, 'lulc_utm32n.tif'),
            'biophysical_table_path': os.path.join(
                DATA, 'biophysical-table.csv'),
            'aoi_vector_path': os.path.join(DATA, 'aoi.gpkg'),
            'cc_method': 'factors',  # or 'intensity'

            # Calibration tool takes modified type (list)
            'ref_eto_table': os.path.join(DATA, 'ref_eto.csv'),

            # t_ref is the reference air temperature. Single value only.
            't_ref': 20,  # from one of Marti's tests TODO: Also test for list of vals
            't_raster_filepaths': glob.glob(os.path.join(DATA, 'T*.tif')),

            't_rasters_table': os.path.join(DATA, 't_rasters.csv'),

            'uhi_maxs': [20],  # tests have scalar 20
            # TODO also test for list of floats, also list of strings.

            # nonstandard args, for the calibration tool
            'station_t_filepath': os.path.join(DATA, 'station-t.csv'),
            'station_locations': os.path.join(DATA, 'station-locations.csv'),
            'station_t_one_day_filepath': os.path.join(DATA, 'station-t-one-day.csv'),
            'num_steps': 2,
            'num_update_logs': 2,
            'initial_solution': '',  # initial values from the model used?  TODO
            'metric': 'RMSE',  # allowed metrics: [R2, MAE, RMSE], optional
            'stepsize': 0.3,
            'exclude_zero_kernel_dist': True,
            'extra_ucm_args': None,  # test key:value parsing, also dict of opts
        }

        with open(self.base_kwargs['ref_eto_table'], 'w') as ref_eto:
            ref_eto.write(textwrap.dedent(
                f"""\
                eto_path,
                {os.path.join(DATA, 'ref_et0.tif')},
                {os.path.join(DATA, 'ref_et1.tif')},
                """))

        with open(self.base_kwargs['t_rasters_table'], 'w') as t_rasters:
            t_rasters.write(textwrap.dedent(
                f"""\
                t_raster_date,t_raster_path,
                23-07-2020,{os.path.join(DATA, '_T0.tif')},
                01-01-2021,{os.path.join(DATA, '_T1.tif')},
                """))

        # Read in the CSV (which is wgs84).  GDAL _can_ read this

        # src_stations has structure {station_name: {x: value, y: value}}
        src_stations = pandas.read_csv(
            os.path.join(DATA, 'station-locations.csv'),
            index_col=0).to_dict(orient='index')

        # station_to_data has structure: {station_name: {date: temp_value}}
        station_t_df = pandas.read_csv(
            os.path.join(DATA, 'station-t.csv'), index_col=0)
        station_t_data = station_t_df.to_dict()

        self.wgs84_srs = osr.SpatialReference()
        self.wgs84_srs.ImportFromEPSG(4326)

        geoms = []
        attributes = []
        fields = {'name': ogr.OFTString}
        for fieldname in station_t_df.index:
            fields[str(fieldname)] = ogr.OFTReal

        for station_name, geom in src_stations.items():
            geoms.append(Point(geom['x'], geom['y']))
            field_data = {
                'name': station_name,
            }

            for date, temp in station_t_data[station_name].items():
                field_data[str(date)] = temp

            attributes.append(field_data)

        self.stations_vector = os.path.join(self.workspace, 'stations.geojson')
        if os.path.exists(self.stations_vector):
            os.remove(self.stations_vector)  # driver won't overwrite files.
        pygeoprocessing.shapely_geometry_to_vector(
            geoms, self.stations_vector, self.wgs84_srs.ExportToWkt(), 'GeoJSON', fields,
            attributes, ogr.wkbPoint)

        invest_args = datastack.extract_parameter_set(
            os.path.join(INVEST_DATA, 'UrbanCoolingModel',
                         'urban_cooling_model_datastack.invest.json')).args

        ref_eto_table_path = os.path.join(
            self.workspace, 'ref_eto_table.csv')
        with open(ref_eto_table_path, 'w') as ref_eto:
            ref_eto.write(textwrap.dedent(
                f"""\
                eto_path,
                {invest_args['ref_eto_raster_path']},
                """))

        # Construct a temperature raster by scaling the LULC's values to the
        # range (18, 24) and then convolving to smooth it out.
        lulc_array = pygeoprocessing.raster_to_numpy_array(
            invest_args['lulc_raster_path']).astype(numpy.float32)
        raster_info = pygeoprocessing.get_raster_info(
            invest_args['lulc_raster_path'])
        nodata_value = raster_info['nodata'][0]
        gt = raster_info['geotransform']
        lulc_min = numpy.min(lulc_array)
        lulc_max = numpy.max(lulc_array)
        rescaled = (
            ((lulc_array - lulc_min) / (lulc_max - lulc_min)) *
            (24-18) + 18)
        convolved = scipy.ndimage.gaussian_filter(rescaled, 3)
        nodata_mask = pygeoprocessing.array_equals_nodata(
            convolved, nodata_value)
        convolved[nodata_mask] = nodata_value
        target_temp_raster = os.path.join(self.workspace, 'temp.tif')
        pygeoprocessing.numpy_array_to_raster(
            convolved, nodata_value, pixel_size=raster_info['pixel_size'],
            origin=(gt[0], gt[3]),
            projection_wkt=raster_info['projection_wkt'],
            target_path=target_temp_raster)

        self.t_rasters_table_path = os.path.join(
            self.workspace, 't_rasters_table.csv')
        with open(self.t_rasters_table_path, 'w') as t_rasters_table:
            t_rasters_table.write(textwrap.dedent(
                f"""\
                t_raster_date,t_raster_path,
                03-11-2025,{target_temp_raster},
                """))

        self.workspace_dir = os.path.join(self.workspace, 'workspace')
        self.calibration_args = {
            'workspace_dir': self.workspace_dir,
            'lulc_raster_path': invest_args['lulc_raster_path'],
            'biophysical_table_path': invest_args['biophysical_table_path'],
            'cc_method': invest_args['cc_method'],
            'ref_eto_table': ref_eto_table_path,
            'aoi_vector_path': invest_args['aoi_vector_path'],
            't_refs': invest_args['t_ref'],
            'uhi_maxs': invest_args['uhi_max'],
            't_rasters_table': self.t_rasters_table_path,
            # skipping t stations vector for now.
            # skipping initial solution so calibrator uses defaults
            # skipping metric so it uses defaults
            # skipping stepsize so it uses defaults
            # Skipping exclude_zero_kernel_dist --> defaults
            # skipping num_steps: defaults
            # skipping num_update_logs: defaults
            # skipping extra_ucm_args: defaults
        }
        if not os.path.exists(self.calibration_args['workspace_dir']):
            os.makedirs(self.calibration_args['workspace_dir'])

    def tearDown(self):
        shutil.rmtree(self.workspace)

    def test_execute(self):
        args = {
            'workspace_dir': self.workspace_dir,
            'lulc_raster_path': self.base_kwargs['lulc_raster_path'],
            'biophysical_table_path': self.base_kwargs['biophysical_table_path'],
            'ref_eto_table': self.base_kwargs['ref_eto_table'],
            'cc_method': 'intensity',
            'aoi_vector_filepath': '',
            't_stations': self.stations_vector,
        }
        plugin.execute(args)

    #def test_execute_with_invest_data(self, t_stations, initial_solution,
    #                                  metric, stepsize,
    #                                  exclude_zero_kernel_dist, num_steps,
    #                                  num_update_logs):
    def test_execute_with_invest_data_and_all_defaults(self):
        plugin.execute(self.calibration_args)
        target_local_dir = os.path.join(CWD, 'test-invest-output')
        if os.path.exists(target_local_dir):
            shutil.rmtree(target_local_dir)
        shutil.copytree(self.workspace, target_local_dir)

    def test_validate(self):
        args = self.base_kwargs.copy()
        args['workspace_dir'] = self.workspace_dir
        self.assertEqual([], plugin.validate(args))

    #@given(n=st.sampled_from(['foo', 'bar', 'baz']),
    #       m=st.sampled_from([1, 2, 3]))
    #def test_hypothesis(self, n, m):
    #    print(n, m)

    def test_vector_to_csvs(self):
        t_stations_vector_path = self.stations_vector
        target_stations_csv = os.path.join(
            self.workspace, 'target_stations.csv')
        target_temps_csv = os.path.join(
            self.workspace, 'target_temps.csv')

        # run the function
        plugin._t_stations_vector_to_csv(
            t_stations_vector_path, target_stations_csv, target_temps_csv)

        for reference_file, new_file in [
                (self.base_kwargs['station_t_filepath'], target_temps_csv),
                (self.base_kwargs['station_locations'], target_stations_csv)]:
            ref_df = pandas.read_csv(reference_file)
            new_df = pandas.read_csv(new_file)
            pandas.testing.assert_frame_equal(ref_df, new_df)


class ParametrizedUCMTests(ParametrizedTestCase):
    """Subclass to allow pytest parametrization."""
    def setUp(self):
        UCMCalibrationPluginTests.setUp(self)

    def tearDown(self):
        UCMCalibrationPluginTests.tearDown(self)

    @parametrize('temps_mode', ['vector', 'rasters'])
    def test_execute_parametrized(self, temps_mode):
        if temps_mode == 'vector':
            stations_vector = os.path.join(
                self.workspace_dir, 'randomized-stations.geojson')
            lulc_raster_info = pygeoprocessing.get_raster_info(
                self.calibration_args['lulc_raster_path'])
            geoms = _randomize_pts_in_bbox(
                lulc_raster_info['bounding_box'], 3, seed=12345)
            dates = ['01-22-2025', '06-20-2025']
            random.seed(100)

            temps = []
            for index, _ in enumerate(geoms):
                values = {'name': f'name_{index}'}
                for date in dates:
                    values[date] = random.uniform(25, 32)
                temps.append(values)

            fields = {'name': ogr.OFTString}
            for date in dates:
                fields[date] = ogr.OFTReal
            pygeoprocessing.shapely_geometry_to_vector(
                geoms, stations_vector, lulc_raster_info['projection_wkt'],
                'GeoJSON', fields, temps, ogr.wkbPoint)

            self.calibration_args['t_stations'] = stations_vector
        elif temps_mode == 'rasters':
            self.calibration_args['t_rasters_table'] = (
                self.t_rasters_table_path)
        plugin.execute(self.calibration_args)
