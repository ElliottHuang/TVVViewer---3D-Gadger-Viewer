import os
from enum import Enum
from trame.app import get_server
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import vtk, vuetify

from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersCore import vtkContourFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonCore import vtkLookupTable
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor, vtkCubeAxesActor, vtkAxesActor
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader
from vtkmodules.vtkFiltersGeneral import vtkWarpVector
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkTextProperty,
)
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkInteractionWidgets import vtkOrientationMarkerWidget
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleSwitch
import vtkmodules.vtkRenderingOpenGL2
import numpy as np

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

DEFAULT_SCALE_FACTOR = 1
ACTIVE_CUBE_AXES_ACTOR = 0
CURRENT_DIRECTORY = os.path.abspath(os.path.dirname(__file__))

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server()
server.client_type = "vue2"
args = server.cli.parse_args()
def get_sub_dir_file(root_dir):
    if not root_dir:
        raise Exception("Base dir should not be empty")
    sub_dir_files = {}
    sub_dir_list = []
    for item in os.listdir(root_dir):
        sub_path = os.path.join(root_dir, item)
        if os.path.isdir(sub_path):
            sub_dir_list.append(item)
            file_list = []
            for f in os.listdir(sub_path):
                f_path = os.path.join(sub_path, f)
                if os.path.isfile(f_path):
                    file_list.append(f)
                elif os.path.isdir(f_path):
                    print("there is dir in ", sub_path)
            sub_dir_files[item] = file_list
        elif os.path.isfile(sub_path):
            print("there is file in ", root_dir)        
    if not sub_dir_list:
        raise Exception("no sub dir")
    elif not sub_dir_files:
        raise Exception("no file in sub dir")
    return sub_dir_list, sub_dir_files

# root_dir = args.dir
root_dir = "./data"
sub_dir_list, sub_dir_files = get_sub_dir_file(root_dir)
cur_sub_dir_idx = 0
cur_file_idx = 0

def getVtkFileName(sub_dir_idx=0, file_idx=0, id=-1):
    if id != -1:
        for key, val in vtk_file_dict.items():
            if val['id'] == id:
                return key
    return sub_dir_files[sub_dir_list[sub_dir_idx]][file_idx]

def getVtkFilePath(sub_dir_idx, file_idx):
    return os.path.join(root_dir, sub_dir_list[sub_dir_idx], getVtkFileName(sub_dir_idx, file_idx))

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

class Representation(Enum):
    Points = 0
    Wireframe = 1
    Surface = 2
    SurfaceWithEdges = 3

class ColorLookupTable(Enum):
    Rainbow = 0
    Inverted_Rainbow = 1
    Greyscale = 2
    Inverted_Greyscale = 3

class Orientation(Enum):
    pos_z = 0
    pos_y = 1
    pos_x = 2
    neg_z = 3
    neg_y = 4
    neg_x = 5

class BG_Color(Enum):
    gray = 0
    black = 1

# -----------------------------------------------------------------------------
# VTK pipeline
# -----------------------------------------------------------------------------

colors = vtkNamedColors()
renderer = vtkRenderer()
renderer.SetBackground(colors.GetColor3d('Black')) # css3 color name
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
# cube_axes = vtkCubeAxesActor()

def warpVectorFilter(reader):
    warpVector = vtkWarpVector()
    warpVector.SetInputConnection(reader.GetOutputPort())
    warpVector.Update()
    return warpVector

# Cube Axes
def CubeAxesActor(mesh_actor):
    cube_axes = vtkCubeAxesActor()
    renderer.AddActor(cube_axes)
    cube_axes.SetBounds(mesh_actor.GetBounds())
    cube_axes.SetCamera(renderer.GetActiveCamera())
    cube_axes.SetXLabelFormat("%6.1f")
    cube_axes.SetYLabelFormat("%6.1f")
    cube_axes.SetZLabelFormat("%6.1f")
    cube_axes.SetFlyModeToOuterEdges()
    renderer.ResetCamera()

def OrientationMarker():
    axes = vtkAxesActor()
    axes.SetAxisLabels(True)
    axes.AxisLabelsOn()
    widget = vtkOrientationMarkerWidget()
    rgba = [0] * 4
    # colors.GetColor('Carrot', rgba)
    # widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(renderWindowInteractor)
    widget.SetViewport(0.0, 0.0, 0.5, 0.5)
    widget.SetEnabled(1)
    widget.InteractiveOff()
    return widget


def calculate_max_length(reader):
    bounds = reader.GetOutput().GetBounds()
    max_vector_length = 0.0
    num_vector = reader.GetNumberOfVectorsInFile()
    if num_vector != 0:
        for i in range(num_vector):
            vector_data = reader.GetOutput().GetPointData().GetArray(0)
            vector_array = np.array(vector_data)
            for v in vector_array:
                vector_magnitude = np.linalg.norm(v)
                max_vector_length = max(max_vector_length, vector_magnitude)
        x_length = bounds[1] - bounds[0]
        y_length = bounds[3] - bounds[2]
        z_length = bounds[5] - bounds[4]
        max_length_of_bounds = max(x_length, y_length, z_length)/2
        if max_vector_length != 0:
            max_length_of_warpbyvector = max_length_of_bounds / max_vector_length
        if(max_length_of_warpbyvector < 1):
            max_length_of_warpbyvector = 1
        else:
            max_length_of_warpbyvector = round(max_length_of_warpbyvector)
        return max_length_of_warpbyvector
    else:
        return 0

def append_scalar(vector_name, reader, path):
    add1 = vector_name + "_magnitude"
    add2 = vector_name + "_x"
    add3 = vector_name + "_y"
    add4 = vector_name + "_z"
    
    with open(path, "r") as vtk_file:
        contents = vtk_file.read()
    vector_data = reader.GetOutput().GetPointData().GetArray(0)
    vector_array = np.array(vector_data)
    contents = contents + "\n" + "SCALARS " + add1 + " double" + "\n" + "LOOKUP_TABLE default"
    for v in vector_array:
        vector_magnitude = np.linalg.norm(v)
        rounded_number = round(vector_magnitude, 7)
        contents = contents + "\n" + str(rounded_number)
    contents = contents + "\n" + "SCALARS " + add2 + " double" + "\n" + "LOOKUP_TABLE default"
    for v in vector_array:
        v_num = v[0]
        contents = contents + "\n" + str(v_num)
    contents = contents + "\n" + "SCALARS " + add3 + " double" + "\n" + "LOOKUP_TABLE default"
    for v in vector_array:
        v_num = v[1]
        contents = contents + "\n" + str(v_num)
    contents = contents + "\n" + "SCALARS " + add4 + " double" + "\n" + "LOOKUP_TABLE default"
    for v in vector_array:
        v_num = v[2]
        contents = contents + "\n" + str(v_num)
    with open(path, "w") as vtk_file:
        vtk_file.write(contents)
    

def process_vtk_file(file_name):
    reader = vtkUnstructuredGridReader()
    if not os.path.exists(file_name):
        raise Exception("file not exist")
    reader.SetFileName(file_name)
    reader.ReadAllScalarsOn() # need since getpointdata().getarray() need
    reader.ReadAllVectorsOn() # same as above
    reader.Update()

    # Calculate the maximum length of the bounding box
    max_length = calculate_max_length(reader)

    # set the name of the scalar data to extract
    num_scalar = reader.GetNumberOfScalarsInFile()
    scalar_list = []
    # print(len(scalar_list))
    for i in range(num_scalar):
        scalar_list.append(reader.GetScalarsNameInFile(i))

    # Set the name of the vector data to extract.
    num_vector = reader.GetNumberOfVectorsInFile()
    vector_list = []
    for i in range(num_vector):
        vector_list.append(reader.GetVectorsNameInFile(i))
    if len(vector_list) == 0:
        raise Exception("no vectors in file")
    else:
        for v in vector_list:
            check_if_added = v + "_magnitude"
            if check_if_added not in scalar_list:
                append_scalar(v, reader, file_name)
    
    num_scalar = reader.GetNumberOfScalarsInFile()
    scalar_list = []
    # print(len(scalar_list))
    for i in range(num_scalar):
        scalar_list.append(reader.GetScalarsNameInFile(i))

    # if len(scalar_list) == 0:
        # raise Exception("no scalars in file")
    #     print("warning: no scalars in file")
    warpVector = warpVectorFilter(reader)

    scalar_range = reader.GetOutput().GetScalarRange()
    # Mesh
    
    mesh_mapper = vtkDataSetMapper()
    mesh_mapper.SetInputConnection(warpVector.GetOutputPort())
    # mesh_mapper.SetLookupTable(mesh_lut)
    mesh_mapper.GetLookupTable().SetRange(scalar_range)
    mesh_mapper.SetUseLookupTableScalarRange(True)
    mesh_mapper.SetScalarVisibility(True)
    mesh_actor = vtkActor()
    mesh_actor.SetMapper(mesh_mapper)
    renderer.AddActor(mesh_actor)

    # Mesh: Setup default representation to surface
    mesh_actor.GetProperty().SetRepresentationToSurface()
    mesh_actor.GetProperty().SetPointSize(1)
    mesh_actor.GetProperty().EdgeVisibilityOff()

    # Mesh: Apply rainbow color map
    # Create a custom lut. The lut is used for both at the mapper and at the scalar_bar
    mesh_lut = mesh_mapper.GetLookupTable()
    mesh_lut.SetHueRange(0.666, 0.0)
    mesh_lut.SetSaturationRange(1.0, 1.0)
    mesh_lut.SetValueRange(1.0, 1.0)
    mesh_lut.Build()
    # Cube axes
    cube_axes = vtkCubeAxesActor()
    renderer.AddActor(cube_axes)
    cube_axes.SetBounds(mesh_actor.GetBounds())
    cube_axes.SetCamera(renderer.GetActiveCamera())
    cube_axes.SetXLabelFormat("%6.1f")
    cube_axes.SetYLabelFormat("%6.1f")
    cube_axes.SetZLabelFormat("%6.1f")
    cube_axes.SetFlyModeToOuterEdges()
    renderer.ResetCamera()
    return mesh_actor, scalar_list, reader, vector_list, warpVector, cube_axes, max_length

vtk_file_dict = {}
init = True
tmp_id = 0
for sub_dir_name, file_list in sub_dir_files.items():
    for f in file_list:
        file_name = os.path.join(root_dir, sub_dir_name, f)
        actor, scalar_list, reader, vector_list, warp_vector, cube_axes, max_length = process_vtk_file(file_name)
        vtk_representation = Representation.Surface.value
        color_by = 0
        color_map = ColorLookupTable.Rainbow.value
        if init: 
            actor.SetVisibility(1)
            cube_axes.SetVisibility(True)
            init = False
        else: 
            actor.SetVisibility(0)
            cube_axes.SetVisibility(False)
        vtk_file_dict[f] = {'id': tmp_id,
                            'actor':actor, 
                            'scalar_list':scalar_list,
                            'reader': reader,
                            'vector_list':vector_list,
                            'warp_vector':warp_vector,
                            'cube_axes':cube_axes,
                            'max_length':max_length,
                            'vtk_representation':vtk_representation,
                            'color_by':color_by,
                            'color_map':color_map
                            }
        tmp_id += 1
        renderer.AddActor(actor)

# create the scalar_bar
scalar_bar = vtkScalarBarActor()
first_vtk = vtk_file_dict[getVtkFileName(0,0)]
lut = first_vtk['actor'].GetMapper().GetLookupTable()
scalar_bar.SetLookupTable(lut)
scalar_bar.GetTitleTextProperty().SetColor(1,0,0)
scalar_bar.DrawTickLabelsOn()
scalar_bar.DrawAnnotationsOn()
scalar_bar.SetNumberOfLabels(10)
if first_vtk['scalar_list'] != []:
    scalar_bar.SetTitle(first_vtk['scalar_list'][0])
else:
    scalar_bar.SetVisibility(0)

renderer.AddActor(scalar_bar)
widget = OrientationMarker()
renderer.ResetCamera()

state, ctrl = server.state, server.controller
state.sub_dir_list = sub_dir_list
state.cur_sub_dir_idx = cur_sub_dir_idx
state.cur_vtk_id = cur_file_idx

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

@state.change("cube_axes_visibility")
def update_cube_axes_visibility(cube_axes_visibility, **kwargs):
    # print(cube_axes_visibility)
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    cur_vtk['cube_axes'].SetVisibility(cube_axes_visibility)
    ctrl.view_update()

@state.change("orientation_marker_visibility")
def update_orientation_marker_visibility(orientation_marker_visibility, **kwargs):
    widget.SetEnabled(orientation_marker_visibility)
    ctrl.view_update()

@state.change("BG_Color")
def update_BG_Color(BG_Color, **kwargs):
    if BG_Color == 0:
        renderer.SetBackground(colors.GetColor3d('Gray'))
    else:
        renderer.SetBackground(colors.GetColor3d('Black'))
    ctrl.view_update()

@state.change("scale_factor")
def update_scale_factor(scale_factor, **kwargs):
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    cur_vtk['warp_vector'].SetScaleFactor(scale_factor)
    ctrl.view_update()

@state.change("mesh_opacity")
def update_mesh_opacity(mesh_opacity, **kwargs):
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    cur_vtk['actor'].GetProperty().SetOpacity(mesh_opacity)
    ctrl.view_update()

@state.change("representation_mode")
def update_representation(representation_mode, **kwargs):    
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    property = cur_vtk['actor'].GetProperty()
    # print(representation_mode)
    if representation_mode == Representation.Points.value:
        property.SetRepresentationToPoints()
        property.SetPointSize(5)
        property.EdgeVisibilityOff()
    elif representation_mode == Representation.Wireframe.value:
        property.SetRepresentationToWireframe()
        property.SetPointSize(1)
        property.EdgeVisibilityOff()
    elif representation_mode == Representation.Surface.value:
        property.SetRepresentationToSurface()
        property.SetPointSize(1)
        property.EdgeVisibilityOff()
    elif representation_mode == Representation.SurfaceWithEdges.value:
        property.SetRepresentationToSurface()
        property.SetPointSize(1)
        property.EdgeVisibilityOn()
    else:
        return
    cur_vtk['vtk_representation'] = representation_mode
    ctrl.view_update()

@state.change("color_map")
def update_colormap(color_map, **kwargs):    
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    lut = cur_vtk['actor'].GetMapper().GetLookupTable()
    if color_map == ColorLookupTable.Rainbow.value:
        lut.SetHueRange(0.666, 0.0)
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
    elif color_map == ColorLookupTable.Inverted_Rainbow.value:
        lut.SetHueRange(0.0, 0.666)
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
    elif color_map == ColorLookupTable.Greyscale.value:
        lut.SetHueRange(0.0, 0.0)
        lut.SetSaturationRange(0.0, 0.0)
        lut.SetValueRange(0.0, 1.0)
    elif color_map == ColorLookupTable.Inverted_Greyscale.value:
        lut.SetHueRange(0.0, 0.666)
        lut.SetSaturationRange(0.0, 0.0)
        lut.SetValueRange(1.0, 0.0)
    else:
        return
    cur_vtk['color_map'] = color_map
    lut.Build()
    ctrl.view_update()

@state.change("color_array_idx")
def update_color_by_name(color_array_idx, **kwargs):
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    scalar_list = cur_vtk['scalar_list']
    reader = cur_vtk['reader']
    cur_vtk['color_by'] = color_array_idx
    # print(scalar_list)
    # print(reader)
    if len(scalar_list) > 0:
        scalar_bar.SetTitle(scalar_list[color_array_idx])
        reader.SetScalarsName(scalar_list[color_array_idx])
        reader.Update()  # Needed because of GetScalarRange
        scalar_range = reader.GetOutput().GetScalarRange()  
        mapper = cur_vtk['actor'].GetMapper()    
        # mapper.SetScalarRange(scalar_range)  not work
        mapper.GetLookupTable().SetRange(scalar_range)
        mapper.SetScalarVisibility(True)
        mapper.SetUseLookupTableScalarRange(True)
        ctrl.view_update()

@state.change("sub_dir_index")
def update_sub_dir_index(sub_dir_index, **kwargs):
    # check if variable exist (getserver() before create global variable will need to check this)
    global cur_sub_dir_idx
    global cur_file_idx
    # print("sub_dir: ", sub_dir_index)
    # print("file: ", file_index)
    def change_sub_dir_index_actor():
        old_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
        new_vtk = vtk_file_dict[getVtkFileName(sub_dir_index, 0)]
        old_vtk['actor'].SetVisibility(0)            
        new_vtk['actor'].SetVisibility(1)
        old_vtk['cube_axes'].SetVisibility(False)            
        new_vtk['cube_axes'].SetVisibility(False)
        new_vtk['warp_vector'].SetScaleFactor(0)
        state.scale_factor = 0
        state.scale_factor_max = new_vtk['max_length']
        
        if new_vtk['scalar_list'] == []:
            scalar_bar.SetVisibility(0)
        else:
            lut = new_vtk['actor'].GetMapper().GetLookupTable()
            scalar_bar.SetLookupTable(lut)
            scalar_bar.SetTitle(new_vtk['scalar_list'][0])
            scalar_bar.SetVisibility(1)
            state.color_array_idx = new_vtk['color_by']
        
        state.cur_vtk_id = new_vtk['id']
        state.representation_mode = new_vtk['vtk_representation'] 
        state.color_map = new_vtk['color_map']
    change_sub_dir_index_actor()
    cur_sub_dir_idx = sub_dir_index
    cur_file_idx = 0
    state.cur_sub_dir_idx = cur_sub_dir_idx
    state.cur_file_idx = cur_file_idx
    state.file_index = 0
    state.cube_axes_visibility = False
    # reader.SetFileName(getVtkFilePath(cur_sub_dir_idx, cur_file_idx))
    # reader.Update()
    # state.representation_mode = Representation.Surface.value
    # state.color_map = ColorLookupTable.Rainbow.value
    ctrl.view_update()

@state.change("file_index")
def update_sub_dir_index(file_index, **kwargs):
    # check if variable exist (getserver() before create global variable will need to check this)
    global cur_sub_dir_idx
    global cur_file_idx
    # print("sub_dir: ", sub_dir_index)
    # print("file: ", file_index)
    def change_file_index_actor():
        old_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
        new_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, file_index)]
        old_vtk['actor'].SetVisibility(0)            
        new_vtk['actor'].SetVisibility(1)
        old_vtk['cube_axes'].SetVisibility(False)            
        new_vtk['cube_axes'].SetVisibility(False)
        new_vtk['warp_vector'].SetScaleFactor(0)
        state.scale_factor = 0
        state.scale_factor_max = new_vtk['max_length']
        
        if new_vtk['scalar_list'] == []:
            scalar_bar.SetVisibility(0)
        else:
            lut = new_vtk['actor'].GetMapper().GetLookupTable()
            scalar_bar.SetLookupTable(lut)
            scalar_bar.SetTitle(new_vtk['scalar_list'][0])
            scalar_bar.SetVisibility(1)
            state.color_array_idx = new_vtk['color_by']
        
        state.cur_vtk_id = new_vtk['id']
        state.representation_mode = new_vtk['vtk_representation'] 
        state.color_map = new_vtk['color_map']
    change_file_index_actor()
    cur_file_idx = file_index
    state.cur_file_idx = cur_file_idx
    state.cube_axes_visibility = False
    # reader.SetFileName(getVtkFilePath(cur_sub_dir_idx, cur_file_idx))
    # reader.Update()
    # state.representation_mode = Representation.Surface.value
    # state.color_map = ColorLookupTable.Rainbow.value
    ctrl.view_update()

def reset_scale_factor():
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    state.scale_factor = 0
    state.set_scale_factor = 0

@state.change("set_scale_factor")
def update_scale_factor(set_scale_factor, **kwargs):
    if not isinstance(set_scale_factor, str):
        rounded_scale_factor = round(set_scale_factor)
    elif set_scale_factor.isdigit():
        rounded_scale_factor = int(set_scale_factor)
    else: 
        try:
            rounded_scale_factor = round(float(set_scale_factor))
        except ValueError:
            # Handle invalid input (e.g., non-numeric strings)
            state.scale_factor = 0  # Or handle the error as appropriate
            rounded_scale_factor = 0
    if rounded_scale_factor < 0:
        state.scale_factor = 0
    elif rounded_scale_factor > state.scale_factor_max:
        state.scale_factor = state.scale_factor_max
    else:
        state.scale_factor = rounded_scale_factor

# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

def tool_bar_icon():
    # with vuetify.VBtn(icon=True, click="$refs.view.resetCamera()"):
    #     vuetify.VIcon("mdi-crop-free")
    
    vuetify.VCheckbox(
        v_model=("cube_axes_visibility", False),
        on_icon="mdi-cube-outline",
        off_icon="mdi-cube-off-outline",
        classes="mx-1",
        hide_details=True,
        dense=True,
    )
    # vuetify.VCheckbox(
    #     v_model=("orientation_marker_visibility", 1),
    #     on_icon="mdi-axis-arrow",
    #     off_icon="mdi-axis-arrow-lock",
    #     classes="mx-1",
    #     hide_details=True,
    #     dense=True,
    # )
    vuetify.VCheckbox(
        v_model=("BG_Color", 1),
        on_icon="mdi-palette",
        off_icon="mdi-palette-outline",
        classes="mx-1",
        hide_details=True,
        dense=True,
    )
    vuetify.VCheckbox(
        v_model="$vuetify.theme.dark",
        on_icon="mdi-lightbulb-off-outline",
        off_icon="mdi-lightbulb-outline",
        classes="mx-1",
        hide_details=True,
        dense=True,
    )
    # vuetify.VBtn(
    #     "RED: X-AXIS",
    #     color="red",
    # )
    # vuetify.VBtn(
    #     "GREEN: Y-AXIS",
    #     color="green",
    # )
    # vuetify.VBtn(
    #     "BLUE: Z-AXIS",
    #     color="blue",
    # )

def drawer_panels(panel_header):
    with vuetify.VExpansionPanel():
        vuetify.VExpansionPanelHeader(
            panel_header,
            classes="grey lighten-1 my-auto py-1 grey--text text--darken-3",
            # style="user-select: none; cursor: pointer",
            hide_details=True,
            dense=False,
        )
        content = vuetify.VExpansionPanelContent(classes="py-2")
    return content

def warp_panel():
    global cur_sub_dir_idx
    global cur_file_idx
    cur_vtk = vtk_file_dict[getVtkFileName(cur_sub_dir_idx, cur_file_idx)]
    with drawer_panels(panel_header="WarpByVector"):
        with vuetify.VRow(classes="pt-2", dense=True):
            with vuetify.VCol():
                vuetify.VSlider(
                # ScaleFactor
                v_model=("scale_factor", 0),
                min=0,
                max=("scale_factor_max", cur_vtk['max_length']),
                step=0.1,
                classes="mt-1",
                hide_details=True,
                dense=True,
                )            
            with vuetify.VCol(cols="4"):
                vuetify.VTextField(
                # display ScaleFactor
                "ScaleFactor",
                v_model=("scale_factor", 0),
                type="number",
                hint="ScaleFactor",
                persistent_hint=True,
                variant="plain",
                dense=True,
                readonly=True,
                classes="mt-1",
                )
            with vuetify.VRow(classes="pt-2", dense=True):
                vuetify.VTextField(
                # set ScaleFactorMax
                "ScaleFactorMax",
                v_model=("set_scale_factor",0),
                type="number",
                label="set ScaleFactorMax",
                placeholder="type number",
                hide_details=True,
                outlined=True,
                dense=True,
                )
                with vuetify.VBtn(icon=True, click=reset_scale_factor):
                    vuetify.VIcon("mdi-restore")

def representation_panel():
     with drawer_panels(panel_header="Representation"):
        vuetify.VSelect(
            v_model=("representation_mode", Representation.Surface.value),
            items=(
                "representations",
                [{"text": r.name, "value": r.value} for r in Representation],
            ),
            label="Representation",
            hide_details=True,
            dense=True,
            outlined=True,
            classes="pt-1",
        )
        with vuetify.VRow(classes="pt-2", dense=True):
            with vuetify.VCol(cols="6"):
                for i in range(len(vtk_file_dict)):
                    vuetify.VSelect(
                        # Color By
                        v_show=(f"{i}==cur_vtk_id"),
                        label="Color by",
                        v_model=("color_array_idx", 0),
                        items=(f"array_list{i}", 
                            [{"text": s, "value": i} for i, s in enumerate(vtk_file_dict[getVtkFileName(id=i)]['scalar_list'])],
                        ),
                        hide_details=True,
                        dense=True,
                        outlined=True,
                        classes="pt-1",
                    )
            with vuetify.VCol(cols="6"):
                vuetify.VSelect(
                    v_model=("color_map", ColorLookupTable.Rainbow.value),
                    items=(
                        "colormaps",
                        [{"text": color.name, "value": color.value} for color in ColorLookupTable],
                    ),
                    label="ColorMap",
                    hide_details=True,
                    dense=True,
                    outlined=True,
                    classes="pt-1",
                )
        vuetify.VSlider(
            # Opacity
            v_model=("mesh_opacity", 1.0),
            min=0,
            max=1,
            step=0.1,
            label="Opacity",
            classes="mt-1",
            hide_details=True,
            dense=True,
        )

def vtk_file_chooser():
    with vuetify.VSlideGroup(v_model=("sub_dir_index",0), show_arrows=True, mandatory=True, classes="mt-2"):
        with vuetify.VSlideItem(v_for=("dir in sub_dir_list",), key=("dir",), v_slot="{ active, toggle }"):
            vuetify.VBtn(
                "{{ dir }}",
                classes="mx-2 mb-1",
                input_value=("active",),
                active_class="primary",
                rounded=True,
                click="toggle"
            )

    for i in range(len(sub_dir_list)):   
        # print(sub_dir_files[sub_dir_list[i]])         
        s = vuetify.VSelect(
            # FileSelect
            v_if=(f"{i}==cur_sub_dir_idx",),
            v_model=("file_index", 0),
            items=(
                f"files{i}",
                [{"text": file_name, "value": idx} for idx, file_name in enumerate(sub_dir_files[sub_dir_list[i]])],
            ),
            label="vtkFile",
            hide_details=True,
            dense=True,
            outlined=True,
            classes="pt-1",
        )

with SinglePageWithDrawerLayout(server) as layout:
    layout.title.set_text("TVVViewer - 3D Gadger Viewer")
    with layout.root:
        with vuetify.VContainer(
            fluid=True,
            classes="pa-0 fill-height",
        ):
            view = vtk.VtkLocalView(renderWindow)
            ctrl.view_update = view.update
            ctrl.view_reset_camera = view.reset_camera
            ctrl.view_widgets_set = view.set_widgets
            view.set_widgets([widget])  # or at constructor
    with layout.toolbar:
        vuetify.VSpacer()
        vuetify.VDivider(vertical=True, classes="mx-2")
        tool_bar_icon()
    with layout.drawer as drawer:
        # drawer components
        drawer.width = 325    
        with vuetify.VExpansionPanels(
            multiple=True,
        ):
            vtk_file_chooser()
            warp_panel()
            representation_panel()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    server.start()
