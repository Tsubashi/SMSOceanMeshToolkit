
# 1. Standard python modules
import sys, logging

# 2. Third party modules
import numpy as np
import geopandas as gpd
import matplotlib
from pyproj import CRS


# 3. Aquaveo modules
# try:
#    from xms.data_objects.parameters import FilterLocation
# except ImportError:  # pragma no cover - optional import
#    FilterLocation = object
from xms.gdal.utilities import gdal_utils as gu
from xms.tool_core import IoDirection, Tool
from xms.grid.ugrid import UGrid as XmUGrid
from xms.constraint.ugrid_builder import UGridBuilder

from xms.tool.utilities.coverage_conversion import polygons_to_shapefile

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

ARG_INPUT_DEM = 8

ARG_MIN_MESH_SIZE = 9
ARG_MAX_EDGE_LENGTH = 10
# make an arg for maximum element bounds based on depth
ARG_MAX_ELEMENT_SIZE_BY_DEPTH = 11
ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS = 12


ARG_GRADATION_RATE = 13

ARG_SIZING_FUNCTION_1 = 14
ARG_NUM_ELEMENTS_PER_SHORELINE = 15
ARG_MAX_ELEMENT_SIZE_NEARSHORE = 16
ARG_NEARSHORE_TOLERANCE = 17

ARG_SAVE_OFF_MEDIAL_AXIS = 18
ARG_SAVE_FILENAME_MEDIAL_AXIS = 19
ARG_LOAD_IN_MEDIAL_AXIS = 20
ARG_LOAD_FILENAME_MEDIAL_AXIS = 21

ARG_SIZING_FUNCTION_2 = 22
ARG_NUM_ELEMENTS_PER_WAVELENGTH = 23
ARG_PERIOD_OF_WAVE = 24

ARG_SIZING_FUNCTION_3 = 25
ARG_DESIRED_TIMESTEP = 26
ARG_MAX_CFL = 27

ARG_FINAL_SIZING_FUNCTION_RASTER = 28

ARG_CLEAN_MESH = 29
ARG_MIN_BOUNDARY_QUALITY = 30
ARG_MIN_PERCENT_DISCONN_AREA = 31
ARG_MAX_LAPLACE_ITER = 32
ARG_MAX_LAPLACE_MOVT_TOL = 33

ARG_ADV_MESH_GENERATION = 34
# number of meshing iterations
ARG_NUM_MESHING_ITERATIONS = 35
ARG_MESHING_PSEUDO_DT = 36
ARG_MESHING_FORCE_FUNCTION = 37

ARG_OUTPUT_UGRID = 38


class OceanMeshUGridTool(Tool):
    """Tool to create a raster of a feature mesh sizing fucntion."""

    def __init__(self):
        """Initializes the class."""
        super().__init__(name="Ocean Mesh UGrid from coverage")

    def enable_arguments(self, arguments):
        """Called to show/hide arguments, change argument values and add new arguments.

        Args:
            arguments(list): The tool arguments.
        """
        # turn off all the arguments by default
        for arg in arguments:
            arg.hide = True

        # show the selection
        arguments[ARG_TYPE_OF_INPUT].hide = False

        # show the input coverage arguments
        if arguments[ARG_TYPE_OF_INPUT].value == "Map coverage":
            # show the map coverage arguments
            arguments[ARG_INPUT_COVERAGE].hide = False

        if arguments[ARG_TYPE_OF_INPUT].value == "Vector file":
            # show the vector file arguments
            arguments[ARG_INPUT_COVERAGE_SHAPEFILE].hide = False

        # make the DEM always visible
        arguments[ARG_INPUT_DEM].hide = False

        # show the maximum element size by depth arguments
        arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].hide = False
        if arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value:
            arguments[ARG_MAX_ELEMENT_SIZES_BY_DEPTH_BOUNDS].hide = False

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
        # hide them if the user doesn't want to modify them
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
            arguments[ARG_NEARSHORE_TOLERANCE].hide = False

            # if the user wants to save off medial axis file, then show the argument
            arguments[ARG_SAVE_OFF_MEDIAL_AXIS].hide = False
            if arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value:
                arguments[ARG_SAVE_FILENAME_MEDIAL_AXIS].hide = False

            # if the user wants to load in a medial axis file, then show the argument
            arguments[ARG_LOAD_IN_MEDIAL_AXIS].hide = False
            if arguments[ARG_LOAD_IN_MEDIAL_AXIS].value:
                arguments[ARG_LOAD_FILENAME_MEDIAL_AXIS].hide = False
            
            # set the save off medial axis to false as we can only load in OR save off
            # if save off is selected, then load in must be false
            if arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value:
                arguments[ARG_LOAD_IN_MEDIAL_AXIS].value = False
                arguments[ARG_LOAD_IN_MEDIAL_AXIS].hide = True
            
            # if load in medial axis is selected, then save off must be false
            if arguments[ARG_LOAD_IN_MEDIAL_AXIS].value:
                arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value = False
                arguments[ARG_SAVE_OFF_MEDIAL_AXIS].hide = True
            
        
        if arguments[ARG_SIZING_FUNCTION_2].text_value == "Wavelength-to-gridscale":
            arguments[ARG_NUM_ELEMENTS_PER_WAVELENGTH].hide = False
            arguments[ARG_PERIOD_OF_WAVE].hide = False
            arguments[ARG_INPUT_DEM].hide = False

        if arguments[ARG_SIZING_FUNCTION_3].text_value == "CFL Timestep Bounding":
            arguments[ARG_DESIRED_TIMESTEP].hide = False
            arguments[ARG_MAX_CFL].hide = False
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
                name="Input Coverage Type",
                description="Specify inputs from map coverage or vector file",
                value="Map coverage",
                choices=["Map coverage", "Vector file"],
            ),
            self.coverage_argument(
                name="Input polgyons from coverage",
                description="Map coverage with polygons",
                optional=True,
            ),
            self.file_argument(
                name="Input polygons as vector file",
                description="Input coverage as a vector file",
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
            self.float_argument(
                name="Moving smoothing window",
                description="Smoothing window applied to polygons (must be odd number or 0)",
                value=5.0,
                optional=True,
            ),
            self.float_argument(
                name="Minimum area multiplier",
                description="Polygons with area < min_area_multi * min. mesh size**2) are removed",
                value=4.0,
                optional=True,
            ),
            self.coverage_argument(
                name="Domain coverage",
                description="Domain coverage as a polygon",
                optional=True,
            ),
            # invert domain
            self.bool_argument(
                name="Invert domain",
                description="Invert the area to-be-meshed (area meshed becomes land if was ocean and vice versa)",
                value=False,
                optional=True,
            ),
            # for DEM
            self.file_argument(
                name="Input Digital Elevation Model (DEM).",
                description="Input DEM for select mesh sizing functions",
                io_direction=IoDirection.INPUT,
                optional=True,
            ),
            # min mesh size
            self.float_argument(
                name="Minimum mesh size",
                description="Minimum mesh size (ft/m)",
                value=1000.0,
            ),
            # max edge length
            self.float_argument(
                name="Maximum mesh size",
                description="Maximum edge length (ft/m)",
                value=99999.0,
                optional=True,
            ),
            # max element size by depth toggl e
            self.bool_argument(
                name="Maximum element size by depth",
                description="Bound maximum element size by depth range",
                value=False,
                optional=True,
            ),
            # max element size by depth bounds
            self.string_argument(
                name="Maximum element size by depth bounds",
                description="Maximum element size by depth bounds (SIZE_BOUND,ELEV_MIN,ELEV_MAX)",
                value="",
                optional=True,
            ),
            # mesh size gradation rate
            self.float_argument(
                name="Approximate rate of change of mesh sizes",
                description="Rate of change",
                value=0.05,
                optional=True,
            ),
            # below are the mesh sizing functions selections
            self.string_argument(
                name="Mesh Sizing Function #1 (required)",
                description="Mesh sizing function #1 (required)",
                value="Distance",
                choices=["Distance", "Feature-size"],
            ),
            # number of elements per shorline width
            self.float_argument(
                name="Number of elements per vector width",
                description="Number of elements per vector width",
                value=3,
                optional=True,
            ),
            self.float_argument(
                name="Maximum element size nearshore (ft/m)",
                description="Maximum element size nearshore (ft/m)",
                value=99999.0,
                optional=True,
            ),
            self.float_argument(
                name="Nearshore tolerance (ft/m)",
                description="Distance from shoreline (ft/m) used to apply the maximum element size nearshore",
                value=1000.0,
                optional=True,
            ),
            # boolean 
            self.bool_argument(
                name="Save off medial axis",
                description="Save off the medial axis",
                value=False,
                optional=True,
            ),
            # save filename
            self.file_argument(
                name="Filename to save medial axis",
                description="Filename to save medial axis",
                io_direction=IoDirection.OUTPUT,
                optional=True,
                value="medial_axis.shp",
            ),
            # boolean 
            self.bool_argument(
                name="Load in medial axis",
                description="Load in a medial axis",
                value=False,
                optional=True,
            ),
            # load filename
            self.file_argument(
                name="Filename to load medial axis",
                description="Filename to load medial axis",
                io_direction=IoDirection.INPUT,
                optional=True,
                value="medial_axis.shp",
            ),
            self.string_argument(
                name="Mesh Sizing Function #2 (optional)",
                description="Mesh sizing function #2 (optional). If selected, specify a DEM",
                value="None",
                optional=True,
                choices=["None", "Wavelength-to-gridscale"],
            ),
            # wave length to grid scale
            self.float_argument(
                name="Number of elements per wavelength",
                description="Number of elements per wavelength",
                value=100,
                optional=True,
            ),
            # period of wave 12.42 hours by default
            self.float_argument(
                name="Period of wave",
                description="Period of wave (seconds)",
                value=12.42 * 3600.0,  # M2 period in seconds
                optional=True,
            ),
            self.string_argument(
                name="Mesh Sizing Function #3 (optional)",
                description="Mesh sizing function #3 (optional). If selected, specify a DEM",
                value="None",
                optional=True,
                choices=["None", "CFL Timestep Bounding"],
            ),
            # cfl timestep bounding
            self.float_argument(
                name="Desired Simulation Timestep",
                description="Maximum allowable timestep (seconds)",
                value=1.0,
                optional=True,
            ),
            self.float_argument(
                name="Maximum Allowable Courant Number",
                description="Maximum allowable Courant number for given timestep",
                value=0.8,
                optional=True,
            ),
            self.raster_argument(
                name="Filename of final mesh sizing function",
                description="The filename of the final mesh sizing function",
                io_direction=IoDirection.OUTPUT,
                # optional=True,
                value="final_sizing_function",
            ),
            self.bool_argument(
                name="Modify mesh cleaning options",
                description="Modify the mesh cleaning options...",
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
                name="Advanced mesh generation options",
                description="Modify mesh generation options",
                value=False,
                optional=True,
            ),
            self.float_argument(
                name="Number of meshing iterations",
                description="Number of meshing iterations",
                value=50,
                optional=True,
            ),
            self.float_argument(
                name="Pseudo timestep for meshing",
                description="Pseudo timestep for meshing",
                value=0.1,
                optional=True,
            ),
            self.string_argument(
                name="Force function for edges of mesh",
                description="Force function used in mesh generation",
                value="bossen_heckbert",
                optional=True,
                choices=["bossen_heckbert", "persson_strang"],
            ),
            # output mesh file
            self.grid_argument(
                name="Filename of mesh",
                description="The filename of the generated mesh",
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
            self.logger.error("Maximum edge length must be greater than 0")
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
        if arguments[ARG_MAX_ELEMENT_SIZE_NEARSHORE].value <= 0:
            self.logger.error("Maximum element size nearshore must be greater than 0")
            return False
        if arguments[ARG_NEARSHORE_TOLERANCE].value <= 0:
            self.logger.error("Nearshore tolerance must be greater than 0")
            return False
        # check if loading in or saving off medial axis
        if arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value:
            if arguments[ARG_SAVE_FILENAME_MEDIAL_AXIS].value == "":
                self.logger.error("Filename to save medial axis must be provided")
                return False
        # make sure either save off or load in medial axis is selected
        if arguments[ARG_LOAD_IN_MEDIAL_AXIS].value:
            if arguments[ARG_LOAD_FILENAME_MEDIAL_AXIS].value == "":
                self.logger.error("Filename to load medial axis must be provided")
                return False
        # load in and save off can't both be true 
        if arguments[ARG_LOAD_IN_MEDIAL_AXIS].value and arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value:
            self.logger.error("Can't load in and save off medial axis at the same time")
            return False
        # check the second sizing function
        if arguments[ARG_SIZING_FUNCTION_2].value == "Wavelength-to-gridscale":
            if arguments[ARG_PERIOD_OF_WAVE].value <= 0:
                self.logger.error("Period of wave must be greater than 0")
                return False
            # ensure there's a DEM
            if arguments[ARG_INPUT_DEM].value == "":
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
        if arguments[ARG_SIZING_FUNCTION_3].value == "CFL Timestep Bounding":
            if arguments[ARG_INPUT_DEM].value == "":
                self.logger.error(
                    "A DEM must be provided for CFL Timestep Bounding sizing function"
                )
                return False
        # if max element size by depth is selected, then the bounds must be provided
        if arguments[ARG_MAX_ELEMENT_SIZE_BY_DEPTH].value:
            # make sure a DEM is provided
            if arguments[ARG_INPUT_DEM].value == "":
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
        time.sleep(10)

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

        elif arguments[ARG_TYPE_OF_INPUT].value == "Map coverage":
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
        nearshore_tolerance = arguments[ARG_NEARSHORE_TOLERANCE].value

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

        if sr.IsGeographic():
            min_mesh_size /= 111000.0  # Approximate meters per degree
            max_mesh_size /= 111000.0  # Approximate meters per degree
            max_size_nearshore /= 111000.0  # Approximate meters per degree
            nearshore_tolerance /= 111000.0  # Approximate meters per degree

            if enforce_max_by_depth:
                size_bound /= 111000.0

            epsg_code = "EPSG:4326"
        else:
            wkt = sr.ExportToWkt()
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
        if arguments[ARG_INPUT_DEM].value != "":

            self.logger.info("Building the DEM inputs...")
            dem = smsom.DEM(
                dem_fname, 
                #ll_ur=region.bbox,
                minimum_resolution=min_mesh_size
            )

        if arguments[ARG_SIZING_FUNCTION_1].value == "Distance":

            self.logger.info("Using a distance sizing function")
            szfx_1 = smsom.distance_sizing_function(
                grid,
                coastal_geometry,
                rate=rate_of_change,
                max_edge_length=max_mesh_size,
            )
        elif arguments[ARG_SIZING_FUNCTION_1].value == "Feature-size":

            self.logger.info("Using feature size mesh sizing function")

            # if the user is not going to load in a medial axis file
            # calculate the medial axis and save it off
            if arguments[ARG_LOAD_IN_MEDIAL_AXIS].value == False:
                self.logger.info("Building medial axis...")
                szfx_1 = smsom.feature_sizing_function(
                    grid,
                    coastal_geometry,
                    number_of_elements_per_width=number_of_elements_per_shoreline,
                    max_edge_length=max_mesh_size,
                    max_element_size_nearshore=max_size_nearshore,
                    nearshore_tolerance=nearshore_tolerance,
                    save_medial_axis=arguments[ARG_SAVE_OFF_MEDIAL_AXIS].value, # save the medial axis medial_axis.gpkg file
                    medial_axis_file=arguments[ARG_SAVE_FILENAME_MEDIAL_AXIS].value, 
            )
            else:
                self.logger.info("Loading in medial axis...")
                # load in the medial axis file from a location specified by the user
                szfx_1 = smsom.feature_sizing_function(
                    grid,
                    coastal_geometry,
                    number_of_elements_per_width=number_of_elements_per_shoreline,
                    max_edge_length=max_mesh_size,
                    max_element_size_nearshore=max_size_nearshore,
                    nearshore_tolerance=nearshore_tolerance,
                    medial_axis_points=arguments[ARG_LOAD_FILENAME_MEDIAL_AXIS].value, 
                )

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
                )
                szfx_1 = smsom.combine_sizing_functions(
                    [szfx_1, szfx_2], operation="min"
                )

        if enforce_max_by_depth:
            self.logger.info("Enforcing maximum element size by depth...")
            szfx_1 = smsom.enforce_mesh_size_bounds_elevation(
                szfx_1, dem, [[min_mesh_size, size_bound, min_bound, max_bound]]
            )

        if arguments[ARG_SIZING_FUNCTION_3].value == "CFL Timestep Bounding":
            self.logger.info("Enforcing CFL timestep in mesh sizing function...")
            szfx_1 = smsom.enforce_CFL_condition(
                szfx_1, dem, desired_timestep, courant_number=max_cfl
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

        self.set_output_grid(b, arguments[ARG_OUTPUT_UGRID])

        self.set_output_raster_file(out_path, out_file)
