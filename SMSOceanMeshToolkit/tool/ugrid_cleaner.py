# 1. Standard python modules
import sys, logging
from typing import Any

# 2. Third party modules
import numpy as np

# 3. Aquaveo modules
from xms.grid.ugrid import UGrid as XmUGrid

from xms.tool_core import IoDirection, Tool
from xms.constraint.ugrid_builder import UGridBuilder

# 4. Local modules
import SMSOceanMeshToolkit as smsom
from xms.tool_core.argument import Argument

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

ARG_INPUT_UGRID = 0
ARG_INPUT_MIN_QUAL = 1
ARG_MIN_PERCENT_DISCONNECTED_AREA = 2
ARG_OUTPUT_UGRID = 3

class UGridCleaner(Tool):
    """
    Topological defects in the UGrid are iteratively removed. 
    
    A valid UGrid defined as having the following properties:
    
        1. the vertices of each triangle are arranged in counterclockwise order;
        2. conformity (a triangle is not allowed to have a vertex of another triangle in its interior); and
        3. traversability (the number of boundary segments is equal to the number of boundary vertices, which guarantees a unique path along the mesh boundary).
    """

    def __init__(self):
        super(UGridCleaner, self).__init__(name='Topologically clean a UGrid')
        
    def initial_arguments(self):
        """Get the initial arguments for the tool.

       Must override.

        Returns:
            (list): A list of the initial tool arguments. 
        """
        arguments = [self.grid_argument(name='Input UGrid/Mesh', description='UGrid/Mesh to clean', io_direction=IoDirection.INPUT),
                     self.float_argument(name='Minimum permitted boundary element quality', value=0.15, description='Minimum permitted quality of the boundary elements'), 
                     self.float_argument(name='Minimum percent treshold for disconnected areas', value=0.05, description='Minimum percent of total area to retain a disconnected area'),
                     self.grid_argument(name='Cleaned Mesh', description='Topologically cleaned mesh', io_direction=IoDirection.OUTPUT), 
        ]
        return arguments

    def run(self, arguments: list[Argument]):
        """Override to run the tool.

        Args:
            arguments (list): The tool arguments.
        """
        #import time 
        #time.sleep(15)
        
        dirty_ugrid = self.get_input_grid(arguments[ARG_INPUT_UGRID].text_value)
        dirty_ugrid = dirty_ugrid.ugrid # convert constrained to unconstrained
        # get points and cells from the input grid
        points = dirty_ugrid.locations
        tris = []
        cnt = 0
        cs = dirty_ugrid.cellstream
        while cnt < len(cs):
            # cell_type = cs[cnt]
            num_nodes = cs[cnt + 1]
            start = cnt + 2
            end = start + num_nodes
            if num_nodes == 3:
                tris.extend(cs[start:end])
            else:
                raise ValueError("Only triangle elements are supported")
            cnt = end

        cells = np.asarray(tris)
        min_element_qual = arguments[ARG_INPUT_MIN_QUAL].value
        min_percent_disconnected_area = arguments[ARG_MIN_PERCENT_DISCONNECTED_AREA].value
        
        points, cells = smsom.mesh_clean(points[:,:2], cells.reshape((-1,3)), min_element_qual, min_percent_disconnected_area, 0, 0)

        # Pass the points and cells to the UGridBuilder
        b = UGridBuilder()
        b.set_is_2d()
        #b.set_unconstrained()

        cellstream = []
        for tri in cells:  # 5 = XMU_TRIANGLE
            cellstream.extend([5, 3, tri[0], tri[1], tri[2]])

        ugrid = XmUGrid(points, cellstream)
        b.set_ugrid(ugrid)
        b = b.build_grid()
        self.set_output_grid(b, arguments[ARG_OUTPUT_UGRID])
