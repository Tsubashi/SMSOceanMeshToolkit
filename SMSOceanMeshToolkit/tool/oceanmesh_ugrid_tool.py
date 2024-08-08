
# 1. Standard python modules
import sys, logging

# 2. Third party modules
import numpy as np
import geopandas as gpd
from pathlib import Path
from pyproj import CRS
from osgeo import ogr


# 3. Aquaveo modules
# try:
#    from xms.data_objects.parameters import FilterLocation
# except ImportError:  # pragma no cover - optional import
#    FilterLocation = object
from xms.gdal.utilities import gdal_utils as gu
from xms.tool_core import IoDirection, Tool
from xms.grid.ugrid import UGrid as XmUGrid
from xms.constraint.ugrid_builder import UGridBuilder
from xms.gdal.vectors import VectorInput
from xms.tool.utilities.coverage_conversion import polygons_to_shapefile, convert_points_to_coverage

# 4. Local modules
import SMSOceanMeshToolkit as smsom

# logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

# matplotlib.use("tkagg")

ARG_TYPE_OF_INPUT = 0
ARG_INPUT_COVERAGE = 1
ARG_INPUT_COVERAGE_SHAPEFILE = 2

# input coverage arguments
ARG_TO_SMOOTH = 3
ARG_SMOOTHING_WINDOW = 4
ARG_MIN_AREA_MULT = 5

ARG_DOMAIN_COVERAGE = 6
ARG_INVERT_DOMAIN = 7

ARG_MIN_MESH_SIZE = 8
ARG_MAX_EDGE_LENGTH = 9
ARG_MAX_ELEMENT_SIZE_BY_DEPTH = 10
ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS = 11

ARG_GRADATION_RATE = 12

ARG_SIZING_FUNCTION_1 = 13
ARG_NUM_ELEMENTS_PER_SHORELINE = 14
ARG_MAX_ELEMENT_SIZE_NEARSHORE = 15

ARG_MEDIAL_AXIS_DIALOG = 16
ARG_MEDIAL_AXIS_INPUT = 17
ARG_MEDIAL_AXIS_OUTPUT = 18

ARG_SIZING_FUNCTION_2 = 19
ARG_NUM_ELEMENTS_PER_WAVELENGTH = 20
ARG_PERIOD_OF_WAVE = 21

ARG_SIZING_FUNCTION_3 = 22
ARG_DESIRED_TIMESTEP = 23
ARG_MAX_CFL = 24

ARG_INPUT_DEM = 25

ARG_CLEAN_MESH = 26
ARG_MIN_BOUNDARY_QUALITY = 27
ARG_MIN_PERCENT_DISCONN_AREA = 28
ARG_MAX_LAPLACE_ITER = 29
ARG_MAX_LAPLACE_MOVT_TOL = 30

ARG_ADV_MESH_GENERATION = 31
ARG_NUM_MESHING_ITERATIONS = 32
ARG_MESHING_PSEUDO_DT = 33
ARG_MESHING_FORCE_FUNCTION = 34

ARG_FINAL_SIZING_FUNCTION_RASTER = 35
ARG_OUTPUT_UGRID = 36


class DummyVariable:
    def __init__(self):
        self.value = None


class OceanMeshUGridTool(Tool):
    """Tool to create a raster of a feature mesh sizing fucntion."""

    def __init__(self):
        """Initializes the class."""
        super().__init__(name="Ocean Mesh UGrid from coverage")
        
    def _push_to_gui(self, file, argument, vector_wkt):
        """Display the vector file composed of points to the map."""
        # extract the vector file name from the path 
        path = Path(file)
        vector_fname = path.name
        vi = VectorInput(file)
        layer_type = vi.layer_type
        display_wkt = gu.add_vertical_to_wkt(self.default_wkt, self.vertical_datum, self.vertical_units)
        # create a dummy variable to hold the values which will be displayed on the map
        dummy_variable = DummyVariable()
        dummy_variable.value = argument
        if layer_type in {ogr.wkbPoint, ogr.wkbPointM, ogr.wkbPointZM, ogr.wkbPoint25D}:
            cov_geometry = vi.get_point_features()
            new_cov = convert_points_to_coverage(cov_geometry, vector_fname, vector_wkt, display_wkt)
            self.set_output_coverage(new_cov, dummy_variable)

    def enable_arguments(self, arguments):
        """Called to show/hide arguments, change argument values and add new arguments.

        Args:
            arguments(list): The tool arguments.
        """
        # turn off all the arguments by default
        for arg in arguments:
            arg.hide = True

        NEED_DEM = False
        
        # show the selection
        arguments[ARG_TYPE_OF_INPUT].hide = False

        # show the input coverage arguments
        if arguments[ARG_TYPE_OF_INPUT].value == "Shoreline coverage":
            arguments[ARG_INPUT_COVERAGE].hide = False

        if arguments[ARG_TYPE_OF_INPUT].value == "Vector file":
            arguments[ARG_INPUT_COVERAGE_SHAPEFILE].hide = False
       
        arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].hide = False
        if arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value:
            arguments[ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS].hide = False
            NEED_DEM = True
            
        # enable the basic arguments
        arguments[ARG_ADV_MESH_GENERATION].hide = False
        if arguments[ARG_ADV_MESH_GENERATION].value:
            arguments[ARG_NUM_MESHING_ITERATIONS].hide = False
            arguments[ARG_MESHING_PSEUDO_DT].hide = False
            arguments[ARG_MESHING_FORCE_FUNCTION].hide = False

        # enable the mesh cleaning arguments
        arguments[ARG_CLEAN_MESH].hide = False
        if arguments[ARG_CLEAN_MESH].value:
            arguments[ARG_MIN_BOUNDARY_QUALITY].hide = False
            arguments[ARG_MIN_PERCENT_DISCONN_AREA].hide = False
            arguments[ARG_MAX_LAPLACE_ITER].hide = False
            arguments[ARG_MAX_LAPLACE_MOVT_TOL].hide = False

        arguments[ARG_TO_SMOOTH].hide = False
        if arguments[ARG_TO_SMOOTH].value:
            arguments[ARG_SMOOTHING_WINDOW].hide = False
            arguments[ARG_MIN_AREA_MULT].hide = False

        # options that always should appears
        arguments[ARG_DOMAIN_COVERAGE].hide = False
        arguments[ARG_INVERT_DOMAIN].hide = False
        arguments[ARG_MIN_MESH_SIZE].hide = False
        arguments[ARG_MAX_EDGE_LENGTH].hide = False
        arguments[ARG_GRADATION_RATE].hide = False

        # enable the sizing function arguments (always should appear)
        arguments[ARG_SIZING_FUNCTION_1].hide = False
        arguments[ARG_SIZING_FUNCTION_2].hide = False
        arguments[ARG_SIZING_FUNCTION_3].hide = False

        # the raster file for the final sizing function
        arguments[ARG_FINAL_SIZING_FUNCTION_RASTER].hide = False

        # the output mesh file
        arguments[ARG_OUTPUT_UGRID].hide = False

        if arguments[ARG_SIZING_FUNCTION_1].text_value == "Feature-size":

            arguments[ARG_NUM_ELEMENTS_PER_SHORELINE].hide = False
            arguments[ARG_MAX_ELEMENT_SIZE_NEARSHORE].hide = False

            # show the option to save off the medial axis
            arguments[ARG_MEDIAL_AXIS_DIALOG].hide = False
            # if the user wants to load in a medial axis file
            if arguments[ARG_MEDIAL_AXIS_DIALOG].value == "Use existing medial axis":

                arguments[ARG_MEDIAL_AXIS_INPUT].hide = False

                arguments[ARG_MEDIAL_AXIS_OUTPUT].value = "None"

            elif arguments[ARG_MEDIAL_AXIS_DIALOG].value == "Compute & output medial axis":
                # if the user wants tto compute and save off the medial axis
                arguments[ARG_MEDIAL_AXIS_OUTPUT].hide = False
                # set what's in ARG_MEDIAL_AXIS_INPUT to None
                arguments[ARG_MEDIAL_AXIS_INPUT].value = "None"
            
        # Toggle for the second sizing function
        if arguments[ARG_SIZING_FUNCTION_2].value:
            arguments[ARG_NUM_ELEMENTS_PER_WAVELENGTH].hide = False
            arguments[ARG_PERIOD_OF_WAVE].hide = False
            NEED_DEM = True

        # Toggle for the third sizing function
        if arguments[ARG_SIZING_FUNCTION_3].value:
            arguments[ARG_DESIRED_TIMESTEP].hide = False
            arguments[ARG_MAX_CFL].hide = False
            NEED_DEM = True

        if NEED_DEM:
            arguments[ARG_INPUT_DEM].hide = False

        return arguments

    def initial_arguments(self):
        """Get initial arguments for tool.

        Must override.

        Returns:
            (list): A list of the initial tool arguments.
        """
        arguments = [
            self.string_argument(
                name="Shoreline/Land polygon format",
                description="Shoreline/Land polygon format",
                value="Shoreline coverage",
                choices=["Shoreline coverage", "Vector file"],
            ),
            self.coverage_argument(
                name="Land polygon coverage",
                description="Land polygon coverage",
                optional=True,
            ),
            self.file_argument(
                name="Land polygon vector file",
                description="Land polygon vector file",
                io_direction=IoDirection.INPUT,
                optional=True,
                value="",
            ),
            self.bool_argument(
                name="Modify options to simplify polygons",
                description="Modify options to simplify polygons",
                value=False,
                optional=True,
            ),
            self.integer_argument(
                name="Moving smoothing window",
                description="Smoothing window applied to polygons (must be odd integer or 0)",
                value=5,
                optional=True,
                min_value=0,
            ),
            self.float_argument(
                name="Minimum area multiplier for polygon removal",
                description="Minimum area multiplier for polygon removal",
                value=4.0,
                optional=True,
                min_value=0.0, 
            ),
            self.coverage_argument(
                name="Domain coverage",
                description="Domain coverage",
                optional=True,
            ),
            self.bool_argument(
                name="Invert domain coverage",
                description="Invert the area to-be-meshed (area meshed becomes land if was water and vice versa)",
                value=False,
                optional=True,
            ),
            self.float_argument(
                name="Minimum mesh size",
                description="Minimum mesh size (ft/m)",
                min_value=0.0, 
            ),
            self.float_argument(
                name="Maximum mesh size",
                description="Maximum mesh size (ft/m)",
                optional=False,
                min_value=0.0,
            ),
            # max element size by depth toggle
            self.bool_argument(
                name="Maximum mesh size by elevation",
                description="Bound maximum mesh size by elevation range",
                value=False,
                optional=True,
            ),
            # max element size by depth bounds
            self.string_argument(
                name="Maximum mesh size by elevation bounds",
                description="Maximum mesh size by elevation bounds (SIZE_BOUND,ELEV_MIN,ELEV_MAX)",
                value="",
                optional=True,
            ),
            # mesh size gradation rate
            self.float_argument(
                name="Approximate rate of change of mesh sizes",
                description="Rate of change",
                value=0.05,
                optional=True,
                min_value=1e-9, 
            ),
            # below are the mesh sizing functions selections
            self.string_argument(
                name="Principal mesh sizing function",
                description="Principal mesh sizing function",
                value="Distance",
                choices=["Distance", "Feature-size"],
            ),
            # number of elements per shorline width
            self.integer_argument(
                name="Number of elements per estimated shoreline width",
                description="Number of elements per estimated shoreline width",
                value=3,
                optional=True,
                min_value=1,
            ),
            self.float_argument(
                name="Maximum mesh size nearshore (ft/m)",
                description="Maximum mesh size nearshore (ft/m)",
                optional=False,
                min_value=0.0,
            ),
            # dialog drop down for medial axis 
            self.string_argument(
                name="Medial axis option",
                description="Medial axis option",
                value="Compute & output medial axis",
                choices=["Use existing medial axis", "Compute & output medial axis"],
            ),
            # if the user wants to load in a medial axis file
            self.file_argument(
                name="Load medial axis file",
                description="Load medial axis file",
                value="",
                io_direction=IoDirection.INPUT,
                optional=True,
            ),
            # if the user wants to save off the medial axis
            self.string_argument(
                name="Save medial axis file",
                description="Save medial axis file",
                value='None',
                optional=True,
            ),
            self.bool_argument(
                name="Use wavelength to gridscale",
                description="Use wavelength to gridscale sizing function",
                value=False,
                optional=True,
            ),
            # wave length to grid scale
            self.integer_argument(
                name="Number of elements per wavelength",
                description="Number of elements per wavelength",
                value=300,
                optional=True,
                min_value=1,
            ),
            # period of wave 12.42 hours by default
            self.float_argument(
                name="Period of wave",
                description="Period of wave (seconds)",
                value=12.42 * 3600.0,  # M2 period in seconds
                optional=True,
                min_value=1e-9,
            ),
            self.bool_argument(
                name="Use CFL timestep bounding",
                description="Use CFL timestep bounding",
                value=False,
                optional=True,
            ),
            self.float_argument(
                name="Desired Simulation Timestep",
                description="Maximum allowable timestep (seconds)",
                value=1.0,
                optional=True,
                min_value=1e-9,
            ),
            self.float_argument(
                name="Maximum Allowable Courant Number",
                description="Maximum allowable Courant number for given timestep",
                value=0.8,
                optional=True,
                min_value=1e-9,
            ),
            # for DEM
            self.file_argument(
                name="Input DEM",
                description="Input DEM for select mesh sizing functions",
                io_direction=IoDirection.INPUT,
                optional=True,
            ),
            self.bool_argument(
                name="Modify mesh cleaning options (advanced)",
                description="Modify the mesh cleaning options (advanced)",
                value=False,
                optional=True,
            ),
            # cleaning arguments go here
            self.float_argument(
                name="Minimum boundary quality (0. to 1.)",
                description="Minimum boundary quality (0. to turn off)",
                value=0.05,
                optional=True,
            ),
            self.float_argument(
                name="Minimum percent disconnected area (0. to 0.50)",
                description="Disconnnected areas with this percent of total area are removed (0. to turn off)",
                value=0.05,
                optional=True,
            ),
            self.float_argument(
                name="Maximum Laplace Iterations",
                description="Maximum number of iterations for Laplace smoothing (0 to turn off)",
                value=20.0,
                optional=True,
            ),
            self.float_argument(
                name="Maximum Laplace Movement Tolerance",
                description="Maximum movement tolerance for Laplace smoothing (0. to turn off)",
                value=0.01,
                optional=True,
            ),
            # boolean for advanced mesh generation
            self.bool_argument(
                name="Mesh generation options (advanced)",
                description="Modify mesh generation options (advanced)",
                value=False,
                optional=True,
            ),
            self.integer_argument(
                name="Number of meshing iterations",
                description="Number of meshing iterations",
                value=50,
                optional=True,
                min_value=2,
            ),
            self.float_argument(
                name="Pseudo timestep for meshing",
                description="Pseudo timestep for meshing",
                value=0.1,
                optional=True,
                min_value=1e-9,
            ),
            self.string_argument(
                name="Force function for edges of mesh",
                description="Force function used in mesh generation",
                value="bossen_heckbert",
                optional=True,
                choices=["bossen_heckbert", "persson_strang"],
            ),
            self.raster_argument(
                name="Output mesh sizing function",
                description="Output mesh sizing function",
                io_direction=IoDirection.OUTPUT,
            ),
            self.grid_argument(
                name="Output mesh",
                description="Output mesh",
                io_direction=IoDirection.OUTPUT,
            ),
        ]

        # depending on the mesh sizing function, enable the arguments
        self.enable_arguments(arguments)

        return arguments

    def _validate_arguments(self, arguments):
        """Called to determine if arguments are valid.

        Args:
            arguments (list): The tool arguments.

        Returns:
            bool: True if the arguments are valid, False if they are not
        """
        if arguments[ARG_MIN_MESH_SIZE].value <= 0:
            self.logger.error("Minimum mesh size must be greater than 0")
            return False

        if arguments[ARG_MAX_EDGE_LENGTH].value <= 0:
            self.logger.error("Maximum mesh size must be greater than 0")
            return False

        if (
            arguments[ARG_GRADATION_RATE].value < 0
            or arguments[ARG_GRADATION_RATE].value > 1
        ):
            self.logger.error("Rate of change must be between 0 and 1")
            return False
        # check if the sizing function options are valid
        if arguments[ARG_SIZING_FUNCTION_1].value not in ["Distance", "Feature-size"]:
            self.logger.error("Invalid sizing function")
            return False
        # check number of elemetns per shoreline
        if arguments[ARG_NUM_ELEMENTS_PER_SHORELINE].value <= 0:
            self.logger.error("Number of elements per shoreline must be greater than 0")
            return False
        # check if None first 
        if arguments[ARG_MAX_ELEMENT_SIZE_NEARSHORE].value != None:
            if arguments[ARG_MAX_ELEMENT_SIZE_NEARSHORE].value <= 0:
                self.logger.error("Maximum element size nearshore must be greater than 0")
                return False
        # check the second sizing function
        if arguments[ARG_SIZING_FUNCTION_2].value:
            if arguments[ARG_PERIOD_OF_WAVE].value <= 0:
                self.logger.error("Period of wave must be greater than 0")
                return False
            # ensure there's a DEM
            if arguments[ARG_INPUT_DEM].value == None:
                self.logger.error(
                    "A DEM must be provided for wavelength-to-gridscale sizing function"
                )
                return False
        if (
            arguments[ARG_SMOOTHING_WINDOW].value % 2 == 0
            or arguments[ARG_SMOOTHING_WINDOW].value < 0
        ):
            self.logger.error("Smoothing window must be an odd number or 0 ")
            return False
        # check min_disconn_area
        if (
            arguments[ARG_MIN_PERCENT_DISCONN_AREA].value < 0
            or arguments[ARG_MIN_PERCENT_DISCONN_AREA].value > 1
        ):
            self.logger.error(
                "Minimum percent disconnected area must be between 0 and 1"
            )
            return False
        # check max_laplace_iter
        if arguments[ARG_MAX_LAPLACE_ITER].value < 0:
            self.logger.error("Maximum Laplace iterations must be greater than 0")
            return False
        # check max_laplace_movt_tol
        if arguments[ARG_MAX_LAPLACE_MOVT_TOL].value < 0:
            self.logger.error(
                "Maximum Laplace movement tolerance must be greater than 0"
            )
            return False
        if arguments[ARG_MAX_CFL].value < 0:
            self.logger.error("Maximum Courant number must be greater than 0")
            return False
        if arguments[ARG_DESIRED_TIMESTEP].value <= 0:
            self.logger.error("Desired timestep must be greater than 0")
            return False
        # if CFL is selected, then the DEM must be provided
        if arguments[ARG_SIZING_FUNCTION_3].value:
            if arguments[ARG_INPUT_DEM].value == None:
                self.logger.error(
                    "A DEM must be provided for CFL Timestep Bounding sizing function"
                )
                return False
        # if max element size by depth is selected, then the bounds must be provided
        if arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value:
            # make sure a DEM is provided
            if arguments[ARG_INPUT_DEM].value == None:
                self.logger.error(
                    "A DEM must be provided for maximum element size by depth"
                )
                return False
        # if the max element size by depth bounds are passed make sure it can be converted to three floats
        if arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value:
            try:
                bounds = arguments[ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS].value.split(
                    ","
                )
                if len(bounds) != 3:
                    self.logger.error(
                        "Invalid number of bounds for maximum element size by depth"
                    )
                    return False
                # check if the bounds can be converted to floats

                for bound in bounds:
                    try:
                        float(bound)
                    except ValueError:
                        self.logger.error(
                            "Invalid bounds for maximum element size by depth"
                        )
                        return False
                bounds = [float(bound) for bound in bounds]
                # ensure the first bound is greater than the min mesh size
                if float(bounds[0]) < arguments[ARG_MIN_MESH_SIZE].value:
                    self.logger.error(
                        f"Bounds ({bound[0]}) for maximum element size by depth must be greater than the minimum mesh size"
                    )
                    return False
                # ensure the second depth bound is less than the third
                if float(bounds[1]) >= float(bounds[2]):
                    self.logger.error(
                        "Min depth bound must be less than the max depth bound"
                    )
                    return False
            except ValueError:
                self.logger.error("Invalid bounds for maximum element size by depth")
                return False

        return True

    def _get_bounding_box(self, polygon_coverage, wkt):
        """
        Helper function to get the bounding box of the domain coverage
        """
        
        tmp_vector = "domain.shp"
        if len(polygon_coverage.polygons) == 0:
            self.logger.error(
                "No polygon found in domain coverage...Please run: Feature Objects --> Build Polygons...aborting"
            )
            return
        polygons_to_shapefile(polygon_coverage.polygons, tmp_vector, wkt)
        c = gpd.read_file(tmp_vector)
        lower_left_x = c.bounds.minx[0]
        lower_left_y = c.bounds.miny[0]
        upper_right_x = c.bounds.maxx[0]
        upper_right_y = c.bounds.maxy[0]

        # make sure coordinates make sense
        if lower_left_x >= upper_right_x:
            self.logger.error("Lower left x must be less than upper right x")
            return

        if lower_left_y >= upper_right_y:
            self.logger.error("Lower left y must be less than upper right y")
            return

        bounding_box = np.array(c.iloc[0].geometry.exterior.coords.xy).T

        return bounding_box

    def run(self, arguments):
        """Override to run the tool.#

        Args:
            arguments (list): The tool arguments.
        """
        import time
        time.sleep(15)

        self._validate_arguments(arguments)

        # determine how the user plans to load the data
        if arguments[ARG_TYPE_OF_INPUT].value == "Vector file":
            COVERAGE_FROM_SHAPEFILE = True
            input_coverage = gpd.read_file(
                arguments[ARG_INPUT_COVERAGE_SHAPEFILE].value
            )

            # check that the input_coverage has a crs
            if input_coverage.crs is None:
                self.logger.error(
                    f"The input coverage specified in the file {arguments[ARG_INPUT_COVERAGE_SHAPEFILE]} must have a CRS. Please add this information and try again. Aborting."
                )
                return

        elif arguments[ARG_TYPE_OF_INPUT].value == "Shoreline coverage":
            COVERAGE_FROM_SHAPEFILE = False
            try:
                input_coverage = self.get_input_coverage(
                    arguments[ARG_INPUT_COVERAGE].value
                )
            except ValueError:
                self.logger.error(
                    "No coverage found...Please run: Feature Objects --> Build Polygons...aborting"
                )
                return

        smooth_bool = arguments[ARG_TO_SMOOTH].value
        smoothing_window = arguments[ARG_SMOOTHING_WINDOW].value
        min_area_mult = arguments[ARG_MIN_AREA_MULT].value

        polygon_coverage = self.get_input_coverage(arguments[ARG_DOMAIN_COVERAGE].value)
        invert_domain = arguments[ARG_INVERT_DOMAIN].value

        dem_fname = arguments[ARG_INPUT_DEM].value

        min_mesh_size = arguments[ARG_MIN_MESH_SIZE].value
        max_mesh_size = arguments[ARG_MAX_EDGE_LENGTH].value

        max_size_nearshore = arguments[ARG_MAX_ELEMENT_SIZE_NEARSHORE].value
        #nearshore_tolerance = arguments[ARG_NEARSHORE_TOLERANCE].value

        rate_of_change = arguments[ARG_GRADATION_RATE].value

        number_of_elements_per_shoreline = arguments[
            ARG_NUM_ELEMENTS_PER_SHORELINE
        ].value
        number_of_elements_per_wavelength = arguments[
            ARG_NUM_ELEMENTS_PER_WAVELENGTH
        ].value
        period_of_wave = arguments[ARG_PERIOD_OF_WAVE].value

        max_cfl = arguments[ARG_MAX_CFL].value
        desired_timestep = arguments[ARG_DESIRED_TIMESTEP].value

        # cleaning arg
        clean_mesh = True  # arguments[ARG_CLEAN_MESH].value
        min_element_qual = arguments[ARG_MIN_BOUNDARY_QUALITY].value
        min_percent_disconn_area = arguments[ARG_MIN_PERCENT_DISCONN_AREA].value
        max_laplace_iter = arguments[ARG_MAX_LAPLACE_ITER].value
        max_laplace_movt_tol = arguments[ARG_MAX_LAPLACE_MOVT_TOL].value

        number_of_meshing_iterations = arguments[ARG_NUM_MESHING_ITERATIONS].value
        psuedo_timestep = arguments[ARG_MESHING_PSEUDO_DT].value
        force_function = arguments[ARG_MESHING_FORCE_FUNCTION].value

        enforce_max_by_depth = arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value
        if enforce_max_by_depth:
            bounds = arguments[ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS].value.split(",")
            # already checked
            min_bound = float(bounds[1])
            max_bound = float(bounds[2])
            size_bound = float(bounds[0])

        sr = None if not self.default_wkt else gu.wkt_to_sr(self.default_wkt)
        if sr is None:
            self.logger.error(
                "This tool requires georeferenced data. A projection must be supplied. Aborting."
            )
            return

        wkt = sr.ExportToWkt()
        
        if sr.IsGeographic():
            min_mesh_size /= 111000.0  # Approximate meters per degree
            max_mesh_size /= 111000.0  # Approximate meters per degree
            # only divide if not None 
            if max_size_nearshore != None:
                max_size_nearshore /= 111000.0  # Approximate meters per degree

            nearshore_tolerance = 1000.0
            nearshore_tolerance /= 111000.0  # Approximate meters per degree

            if enforce_max_by_depth:
                size_bound /= 111000.0
            
            epsg_code = "EPSG:4326"
        else:

            nearshore_tolerance = 1000.0 

            # Create a CRS object from the WKT string
            crs = CRS.from_wkt(wkt)

            # Get the EPSG code
            epsg_code = crs.to_epsg()

        # This method exports the polygons to a shapefile as vector file
        # check if polygon_coverage has polygons
        if polygon_coverage is None:
            self.logger.error(
                "No polygons found in domain coverage...Please run: Feature Objects --> Build Polygons...aborting"
            )
            return

        bounding_box = self._get_bounding_box(polygon_coverage, epsg_code)

        if not COVERAGE_FROM_SHAPEFILE:
            if len(input_coverage.polygons) == 0:
                self.logger.error(
                    "No polygons found in coverage...Please run: Feature Objects --> Build Polygons...aborting"
                )
                return

        if COVERAGE_FROM_SHAPEFILE:
            self.logger.info("Creating mesh from vector file...")
            tmp_vector = arguments[ARG_INPUT_COVERAGE_SHAPEFILE].value
        else:
            self.logger.info("Creating a mesh from a map coverage...")
            tmp_vector = "polys.shp"
            polygons_to_shapefile(input_coverage.polygons, tmp_vector, epsg_code)

        # Oceanmesh specific code below
        region = smsom.Region(bounding_box, crs=epsg_code)
        grid = smsom.Grid(region, dx=min_mesh_size)

        # if smoothing_window is 0, then make it 1 so the smoothing function is disabled
        if smoothing_window == 0:
            smoothing_window = 1

        self.logger.info("Building the coastal geometry...")
        coastal_geometry = smsom.CoastalGeometry(
            tmp_vector,
            bounding_box,
            min_mesh_size,
            crs=epsg_code,
            smooth_shoreline=smooth_bool,
            smoothing_approach="moving_window",
            smoothing_window=int(smoothing_window),
            minimum_area_mult=min_area_mult,
        )
        # determine how many sizing functions are being used
        SZFX_COUNT = 1

        # if the dem is available, then build the first sizing function
        if arguments[ARG_INPUT_DEM].value != None:

            self.logger.info("Building the DEM inputs...")
            dem = smsom.DEM(
                dem_fname, 
                ll_ur=region.bbox,
                minimum_resolution=min_mesh_size, 
                target_crs=epsg_code,
                logger=self.logger,
            )

        if arguments[ARG_SIZING_FUNCTION_1].value == "Distance":

            self.logger.info("Using a distance sizing function")
            szfx_1 = smsom.distance_sizing_function(
                grid,
                coastal_geometry,
                rate=rate_of_change,
                max_edge_length=max_mesh_size,
                logger=self.logger,
            )
        elif arguments[ARG_SIZING_FUNCTION_1].value == "Feature-size":

            self.logger.info("Using feature size mesh sizing function")

            if arguments[ARG_MEDIAL_AXIS_OUTPUT].value != "None":
                # compute the medial axis and save it off
                self.logger.info("Building medial axis...")
                maf = arguments[ARG_MEDIAL_AXIS_OUTPUT].value
                szfx_1 = smsom.feature_sizing_function(
                    grid,
                    coastal_geometry,
                    number_of_elements_per_width=float(number_of_elements_per_shoreline),
                    max_edge_length=max_mesh_size,
                    max_element_size_nearshore=max_size_nearshore,
                    nearshore_tolerance=nearshore_tolerance,
                    save_medial_axis=True, 
                    medial_axis_file=maf, 
                    logger=self.logger,
                )
                # Display the medial axis on the map since it was saved
                self._push_to_gui(arguments[ARG_MEDIAL_AXIS_OUTPUT].value, "approximate_medial_axis", wkt)

            else:
                self.logger.info("Loading in medial axis...")
                maf = arguments[ARG_MEDIAL_AXIS_INPUT].value
                # load in the medial axis file from a location specified by the user
                szfx_1 = smsom.feature_sizing_function(
                    grid,
                    coastal_geometry,
                    number_of_elements_per_width=float(number_of_elements_per_shoreline),
                    max_edge_length=max_mesh_size,
                    max_element_size_nearshore=max_size_nearshore,
                    nearshore_tolerance=nearshore_tolerance,
                    save_medial_axis=False,
                    medial_axis_points=maf, 
                    logger=self.logger,
                )
                # file is being loaded so it must exist, lets display it
                self._push_to_gui(arguments[ARG_MEDIAL_AXIS_INPUT].value, "approximate_medial_axis", wkt)

                

        if arguments[ARG_SIZING_FUNCTION_2].value != "None":
            if arguments[ARG_SIZING_FUNCTION_2].value == "Wavelength-to-gridscale":

                self.logger.info(
                    "Building a wavelength-to-gridscale sizing function..."
                )
                SZFX_COUNT += 1
                szfx_2 = smsom.wavelength_sizing_function(
                    grid,
                    dem,
                    wl=number_of_elements_per_wavelength,
                    max_edge_length=max_mesh_size,
                    period=period_of_wave,
                    logger=self.logger,
                )
                szfx_1 = smsom.combine_sizing_functions(
                    [szfx_1, szfx_2], operation="min"
                )

        if enforce_max_by_depth:
            self.logger.info("Enforcing maximum element size by depth...")
            szfx_1 = smsom.enforce_mesh_size_bounds_elevation(
                szfx_1, dem, [[min_mesh_size, size_bound, min_bound, max_bound]], logger=self.logger
            )

        if arguments[ARG_SIZING_FUNCTION_3].value == "CFL Timestep Bounding":
            self.logger.info("Enforcing CFL timestep in mesh sizing function...")
            szfx_1 = smsom.enforce_CFL_condition(
                szfx_1, dem, desired_timestep, courant_number=max_cfl, logger=self.logger
            )

        self.logger.info("Enforcing mesh size gradation...")
        szfx_1 = smsom.enforce_mesh_gradation(szfx_1, gradation=rate_of_change)

        # write it to a tif file
        out_path = self.get_output_raster(
            arguments[ARG_FINAL_SIZING_FUNCTION_RASTER].value
        )

        # convert the xarray to a raster
        ds = szfx_1.to_xarray()

        if not arguments[ARG_FINAL_SIZING_FUNCTION_RASTER].value.endswith(".tif"):
            out_file = arguments[ARG_FINAL_SIZING_FUNCTION_RASTER].value + ".tif"
        else:
            out_file = arguments[ARG_FINAL_SIZING_FUNCTION_RASTER].value

        ds.rio.to_raster(out_file)

        # Build a signed distance function (SDF)
        self.logger.info("Building the signed distance function...")
        sdf = smsom.signed_distance_function(coastal_geometry, invert=invert_domain)

        # Using the szfx and the sdf build the mesh
        self.logger.info("Generating the mesh...")
        points, cells = smsom.generate_mesh(
            domain=sdf,
            edge_length=szfx_1,
            max_iter=int(number_of_meshing_iterations),
            pseudo_dt=psuedo_timestep,
            force_function=force_function,
            logger=self.logger,
        )

        if clean_mesh:

            self.logger.info("Cleaning the mesh...")
            # clean up the mesh
            points, cells = smsom.mesh_clean(
                points,
                cells,
                min_element_qual=min_element_qual,
                min_percent_disconnected_area=min_percent_disconn_area,
                max_iter=int(max_laplace_iter),
                tol=max_laplace_movt_tol,
            )

        self.logger.info("Finished meshing...")

        # Pass the points and cells to the UGridBuilder
        b = UGridBuilder()
        b.set_is_2d()
        # b.set_unconstrained()

        cellstream = []
        for tri in cells:  # 5 = XMU_TRIANGLE
            cellstream.extend([5, 3, tri[0], tri[1], tri[2]])

        ugrid = XmUGrid(points, cellstream)
        b.set_ugrid(ugrid)
        b = b.build_grid()

        self.set_output_grid(b, arguments[ARG_OUTPUT_UGRID], force_ugrid=False)

        self.set_output_raster_file(out_path, out_file)
