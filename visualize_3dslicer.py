import slicer
import vtk
import sys
import os
import argparse

_default_array_name = "label"
_default_color_map_name = "Labels"
_default_opacity = 1.0
_default_layout_id = 4
_default_range_min = None
_default_range_max = None
_default_reset_view = True
_default_capture_path = None
_default_rot_x = 0.0
_default_rot_y = 0.0
_default_rot_z = 0.0
_default_zoom = 1.0
_default_pan_x = 0.0
_default_pan_y = 0.0
_default_parallel = False
_default_fov = 30.0
_default_capture_size = (1024, 768)

# Module-level offscreen resources (created lazily for headless mode)
_offscreen_rw = None
_offscreen_renderer = None


def _has_gui():
    lm = slicer.app.layoutManager()
    return lm is not None and lm.threeDWidget(0) is not None


def _get_render_window_and_renderer():
    """Return (renderWindow, renderer) from GUI or offscreen."""
    if _has_gui():
        view = slicer.app.layoutManager().threeDWidget(0).threeDView()
        rw = view.renderWindow()
        renderer = rw.GetRenderers().GetFirstRenderer()
        return rw, renderer
    return _offscreen_rw, _offscreen_renderer


def _setup_offscreen(model_node, capture_size=_default_capture_size):
    """Create offscreen render window and populate with model actors."""
    global _offscreen_rw, _offscreen_renderer

    _offscreen_rw = vtk.vtkRenderWindow()
    _offscreen_rw.SetOffScreenRendering(True)
    _offscreen_rw.SetSize(*capture_size)

    _offscreen_renderer = vtk.vtkRenderer()
    _offscreen_renderer.SetBackground(1, 1, 1)
    _offscreen_rw.AddRenderer(_offscreen_renderer)

    dn = model_node.GetDisplayNode()
    polydata = model_node.GetPolyData()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray(dn.GetActiveScalarName())
    mapper.SetScalarRange(dn.GetScalarRange())

    color_node = dn.GetColorNode()
    if color_node:
        lut = color_node.GetScalarsToColors()
        if lut:
            mapper.SetLookupTable(lut)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(dn.GetOpacity())
    _offscreen_renderer.AddActor(actor)
    _offscreen_renderer.ResetCamera()


def load_model(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    node = slicer.util.loadModel(path)
    if node is None:
        raise RuntimeError("Failed to load model")

    if not node.GetDisplayNode():
        node.CreateDefaultDisplayNodes()

    return node


def apply_colormap(
    model_node,
    array_name,
    color_map_name,
    opacity,
    range_min,
    range_max,
):
    polydata = model_node.GetPolyData()
    point_data = polydata.GetPointData()

    if not point_data.HasArray(array_name):
        raise ValueError(f"Array '{array_name}' not found")

    point_data.SetActiveScalars(array_name)

    display_node = model_node.GetDisplayNode()
    display_node.SetActiveScalar(array_name, vtk.vtkAssignAttribute.POINT_DATA)
    display_node.SetScalarVisibility(True)

    color_node = slicer.util.getNode(color_map_name)
    if not color_node:
        raise ValueError(f"Color map '{color_map_name}' not found")

    display_node.SetAndObserveColorNodeID(color_node.GetID())

    data_min, data_max = point_data.GetArray(array_name).GetRange()

    if range_min is None:
        range_min = data_min
    if range_max is None:
        range_max = data_max

    # Disable auto scalar range so manual range is not overwritten
    display_node.SetAutoScalarRange(False)
    display_node.SetScalarRange(range_min, range_max)
    display_node.SetOpacity(opacity)

    # Notify VTK pipeline that display properties have changed
    display_node.Modified()


def setup_layout(layout_id):
    if _has_gui():
        slicer.app.layoutManager().setLayout(layout_id)


def reset_3d_view():
    rw, renderer = _get_render_window_and_renderer()
    if _has_gui():
        slicer.app.layoutManager().threeDWidget(0).threeDView().resetFocalPoint()
    elif renderer:
        renderer.ResetCamera()


def adjust_camera(rot_x, rot_y, rot_z, zoom, pan_x, pan_y, parallel, fov):
    rw, renderer = _get_render_window_and_renderer()
    camera = renderer.GetActiveCamera()

    # Projection mode
    if parallel:
        camera.SetParallelProjection(True)
        # Recompute ParallelScale to fit the scene after switching projection
        renderer.ResetCamera()
    else:
        camera.SetParallelProjection(False)
        camera.SetViewAngle(fov)

    # Rotation
    if rot_x != 0.0:
        camera.Elevation(rot_x)
    if rot_y != 0.0:
        camera.Azimuth(rot_y)
    if rot_z != 0.0:
        camera.Roll(rot_z)

    camera.OrthogonalizeViewUp()

    # Zoom (Dolly for perspective, ParallelScale for parallel)
    if zoom != 1.0:
        if camera.GetParallelProjection():
            camera.SetParallelScale(camera.GetParallelScale() / zoom)
        else:
            camera.Dolly(zoom)

    # Pan (shift focal point in camera-local X/Y)
    if pan_x != 0.0 or pan_y != 0.0:
        fp = list(camera.GetFocalPoint())
        pos = list(camera.GetPosition())
        vup = list(camera.GetViewUp())

        # Camera forward vector
        fwd = [fp[i] - pos[i] for i in range(3)]
        # Camera right vector = forward x up
        right = [
            fwd[1] * vup[2] - fwd[2] * vup[1],
            fwd[2] * vup[0] - fwd[0] * vup[2],
            fwd[0] * vup[1] - fwd[1] * vup[0],
        ]
        # Normalize right
        mag = sum(r * r for r in right) ** 0.5
        if mag > 0:
            right = [r / mag for r in right]

        for i in range(3):
            shift = right[i] * pan_x + vup[i] * pan_y
            fp[i] += shift
            pos[i] += shift

        camera.SetFocalPoint(*fp)
        camera.SetPosition(*pos)

    renderer.ResetCameraClippingRange()


def save_3d_capture(output_path: str):
    if output_path is None:
        return

    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    rw, renderer = _get_render_window_and_renderer()

    if _has_gui():
        slicer.app.processEvents()
        threeDView = slicer.app.layoutManager().threeDWidget(0).threeDView()
        threeDView.forceRender()

    rw.Render()

    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(rw)
    windowToImageFilter.SetInputBufferTypeToRGBA()
    windowToImageFilter.ReadFrontBufferOff()
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output_path)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()

    print(f"Saved capture: {output_path}")


def main(
    vtk_path,
    array_name=None,
    color_map_name=None,
    opacity=None,
    range_min=None,
    range_max=None,
    layout_id=None,
    reset_view=None,
    capture_path=None,
    rot_x=None,
    rot_y=None,
    rot_z=None,
    zoom=None,
    pan_x=None,
    pan_y=None,
    parallel=None,
    fov=None,
):
    if array_name is None:
        array_name = _default_array_name
    if color_map_name is None:
        color_map_name = _default_color_map_name
    if opacity is None:
        opacity = _default_opacity
    if layout_id is None:
        layout_id = _default_layout_id
    if range_min is None:
        range_min = _default_range_min
    if range_max is None:
        range_max = _default_range_max
    if reset_view is None:
        reset_view = _default_reset_view
    if capture_path is None:
        capture_path = _default_capture_path
    if rot_x is None:
        rot_x = _default_rot_x
    if rot_y is None:
        rot_y = _default_rot_y
    if rot_z is None:
        rot_z = _default_rot_z
    if zoom is None:
        zoom = _default_zoom
    if pan_x is None:
        pan_x = _default_pan_x
    if pan_y is None:
        pan_y = _default_pan_y
    if parallel is None:
        parallel = _default_parallel
    if fov is None:
        fov = _default_fov

    setup_layout(layout_id)

    model_node = load_model(vtk_path)

    apply_colormap(
        model_node,
        array_name,
        color_map_name,
        opacity,
        range_min,
        range_max,
    )

    # For headless mode: build offscreen VTK pipeline from applied display settings
    if not _has_gui():
        _setup_offscreen(model_node)

    if reset_view:
        reset_3d_view()

    adjust_camera(rot_x, rot_y, rot_z, zoom, pan_x, pan_y, parallel, fov)

    if _has_gui():
        slicer.app.processEvents()

    if capture_path is not None:
        save_3d_capture(capture_path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Operate 3D Slicer auto visualization tool"
    )

    parser.add_argument("vtk_path")

    parser.add_argument("--array", default=_default_array_name)
    parser.add_argument("--colormap", default=_default_color_map_name)
    parser.add_argument("--opacity", type=float, default=_default_opacity)
    parser.add_argument("--range-min", type=float, default=_default_range_min)
    parser.add_argument("--range-max", type=float, default=_default_range_max)
    parser.add_argument("--layout", type=int, default=_default_layout_id)
    parser.add_argument("--no-reset-view", action="store_true")
    parser.add_argument("--capture", default=_default_capture_path)
    parser.add_argument("--rotx", type=float, default=_default_rot_x,
                        help="Camera elevation (pitch) in degrees")
    parser.add_argument("--roty", type=float, default=_default_rot_y,
                        help="Camera azimuth (yaw) in degrees")
    parser.add_argument("--rotz", type=float, default=_default_rot_z,
                        help="Camera roll in degrees")
    parser.add_argument("--zoom", type=float, default=_default_zoom,
                        help="Zoom factor (1.0=no change, 2.0=2x closer)")
    parser.add_argument("--pan-x", type=float, default=_default_pan_x,
                        help="Pan focal point horizontally (world units)")
    parser.add_argument("--pan-y", type=float, default=_default_pan_y,
                        help="Pan focal point vertically (world units)")
    parser.add_argument("--parallel", action="store_true", default=_default_parallel,
                        help="Use parallel (orthographic) projection")
    parser.add_argument("--fov", type=float, default=_default_fov,
                        help="Field of view in degrees (perspective mode)")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    main(
        vtk_path=args.vtk_path,
        array_name=args.array,
        color_map_name=args.colormap,
        opacity=args.opacity,
        range_min=args.range_min,
        range_max=args.range_max,
        layout_id=args.layout,
        reset_view=not args.no_reset_view,
        capture_path=args.capture,
        rot_x=args.rotx,
        rot_y=args.roty,
        rot_z=args.rotz,
        zoom=args.zoom,
        pan_x=args.pan_x,
        pan_y=args.pan_y,
        parallel=args.parallel,
        fov=args.fov,
    )