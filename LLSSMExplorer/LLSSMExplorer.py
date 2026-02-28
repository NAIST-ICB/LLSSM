"""
LLSSMExplorer - 3D Slicer Scripted Module
Real-time LLSSM (Lower Limb Statistical Shape Model) deformation viewer.

Load an LLSSM .npz file and interactively explore PCA shape modes
with a slider (-3σ to +3σ).
"""

import os
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
    ScriptedLoadableModuleLogic,
)
from slicer.util import VTKObservationMixin
import ctk
import qt


# ============================================================
# Module registration
# ============================================================
class LLSSMExplorer(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "LLSSM Explorer"
        self.parent.categories = ["Shape Analysis"]
        self.parent.dependencies = []
        self.parent.contributors = ["Auto-generated LLSSM Explorer"]
        self.parent.helpText = (
            "Interactive PCA shape mode explorer for LLSSM .npz models.\n"
            "Load a model, choose a component, and drag the beta slider."
        )
        self.parent.acknowledgementText = ""


# ============================================================
# Widget (GUI)
# ============================================================
class LLSSMExplorerWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = LLSSMExplorerLogic()

        # --- Collapsible: I/O ---
        io_section = ctk.ctkCollapsibleButton()
        io_section.text = "Load LLSSM"
        self.layout.addWidget(io_section)
        io_layout = qt.QFormLayout(io_section)

        # File path
        self.path_edit = ctk.ctkPathLineEdit()
        self.path_edit.filters = ctk.ctkPathLineEdit.Files
        self.path_edit.nameFilters = ["NPZ files (*.npz)"]
        self.path_edit.settingKey = "LLSSMExplorer/NpzPath"
        initial_npz = self._resolve_initial_npz_path()
        if initial_npz:
            self.path_edit.currentPath = initial_npz
        io_layout.addRow("NPZ File:", self.path_edit)

        # Load button
        self.load_button = qt.QPushButton("Load LLSSM")
        self.load_button.toolTip = "Load LLSSM from the specified NPZ file"
        io_layout.addRow(self.load_button)

        # Status
        self.status_label = qt.QLabel("No model loaded.")
        io_layout.addRow("Status:", self.status_label)

        # --- Collapsible: Controls ---
        ctrl_section = ctk.ctkCollapsibleButton()
        ctrl_section.text = "Shape Controls"
        self.layout.addWidget(ctrl_section)
        ctrl_layout = qt.QFormLayout(ctrl_section)

        # Component selector
        self.component_selector = qt.QComboBox()
        self.component_selector.enabled = False
        ctrl_layout.addRow("Component:", self.component_selector)

        # Beta range (min / max separate)
        range_box = qt.QHBoxLayout()
        self.range_min_spin = ctk.ctkDoubleSpinBox()
        self.range_min_spin.minimum = -10.0
        self.range_min_spin.maximum = 0.0
        self.range_min_spin.value = -3.0
        self.range_min_spin.singleStep = 0.5
        self.range_min_spin.decimals = 1
        self.range_min_spin.suffix = " σ"
        self.range_min_spin.prefix = ""
        range_box.addWidget(self.range_min_spin)

        range_label = qt.QLabel(" to ")
        range_box.addWidget(range_label)

        self.range_max_spin = ctk.ctkDoubleSpinBox()
        self.range_max_spin.minimum = 0.0
        self.range_max_spin.maximum = 10.0
        self.range_max_spin.value = 3.0
        self.range_max_spin.singleStep = 0.5
        self.range_max_spin.decimals = 1
        self.range_max_spin.suffix = " σ"
        range_box.addWidget(self.range_max_spin)

        ctrl_layout.addRow("Range:", range_box)

        # Beta slider
        self.beta_slider = ctk.ctkSliderWidget()
        self.beta_slider.singleStep = 0.01
        self.beta_slider.minimum = -3.0
        self.beta_slider.maximum = 3.0
        self.beta_slider.value = 0.0
        self.beta_slider.decimals = 2
        self.beta_slider.tracking = True  # update while dragging
        self.beta_slider.enabled = False
        ctrl_layout.addRow("Beta (σ):", self.beta_slider)

        # Active betas summary (read-only text showing all non-zero betas)
        self.betas_display = qt.QTextEdit()
        self.betas_display.readOnly = True
        self.betas_display.setMaximumHeight(60)
        self.betas_display.setStyleSheet("font-family: monospace; font-size: 11px;")
        self.betas_display.plainText = "(no model loaded)"
        ctrl_layout.addRow("Active:", self.betas_display)

        # Reset button
        self.reset_button = qt.QPushButton("Reset to Mean")
        self.reset_button.enabled = False
        ctrl_layout.addRow(self.reset_button)

        # --- Collapsible: Display ---
        disp_section = ctk.ctkCollapsibleButton()
        disp_section.text = "Display"
        disp_section.collapsed = True
        self.layout.addWidget(disp_section)
        disp_layout = qt.QFormLayout(disp_section)

        # Color table selector
        self.color_selector = qt.QComboBox()
        self._color_table_names = [
            "GenericAnatomyColors",
            "GenericColors",
            "MediumChartColors",
            "SPL-BrainAtlas-ColorFile",
            "AbdomenColors",
            "64Color-Nonsemantic",
            "Slicer3_2010_Label_Colors",
            "Random",
            "Warm1", "Warm2", "Warm3",
            "Cool1", "Cool2", "Cool3",
            "WarmShade1", "WarmShade2", "WarmShade3",
            "ColdToHotRainbow",
            "PET-Heat",
            "Grey",
        ]
        for name in self._color_table_names:
            self.color_selector.addItem(name)
        # Resolve initial color table from args / env / default
        initial_color = self._resolve_initial_color_table()
        if initial_color and initial_color in self._color_table_names:
            self.color_selector.setCurrentIndex(
                self._color_table_names.index(initial_color)
            )
        else:
            self.color_selector.setCurrentIndex(0)
        disp_layout.addRow("Color Table:", self.color_selector)

        # Dim opacity for unselected structures
        self.dim_opacity_spin = qt.QSpinBox()
        self.dim_opacity_spin.minimum = 0
        self.dim_opacity_spin.maximum = 255
        self.dim_opacity_spin.value = 40
        self.dim_opacity_spin.suffix = " / 255"
        disp_layout.addRow("Dim opacity:", self.dim_opacity_spin)

        # --- Collapsible: Structures ---
        struct_section = ctk.ctkCollapsibleButton()
        struct_section.text = "Structures"
        struct_section.collapsed = True
        self.layout.addWidget(struct_section)
        struct_layout = qt.QVBoxLayout(struct_section)

        # Quick select row
        quick_box = qt.QHBoxLayout()
        self.select_all_btn = qt.QPushButton("All")
        self.select_none_btn = qt.QPushButton("None")
        quick_box.addWidget(self.select_all_btn)
        quick_box.addWidget(self.select_none_btn)
        struct_layout.addLayout(quick_box)

        # Filter toggles — Side
        side_box = qt.QHBoxLayout()
        side_box.addWidget(qt.QLabel("Side:"))
        self.filter_left = qt.QCheckBox("Left")
        self.filter_right = qt.QCheckBox("Right")
        self.filter_center = qt.QCheckBox("Center")
        self.filter_left.setChecked(True)
        self.filter_right.setChecked(True)
        self.filter_center.setChecked(True)
        side_box.addWidget(self.filter_left)
        side_box.addWidget(self.filter_center)
        side_box.addWidget(self.filter_right)
        struct_layout.addLayout(side_box)

        # Filter toggles — Region
        region_box = qt.QHBoxLayout()
        region_box.addWidget(qt.QLabel("Region:"))
        self.filter_hip_thigh = qt.QCheckBox("Hip/Thigh")
        self.filter_lower_leg = qt.QCheckBox("Lower Leg")
        self.filter_hip_thigh.setChecked(True)
        self.filter_lower_leg.setChecked(True)
        region_box.addWidget(self.filter_hip_thigh)
        region_box.addWidget(self.filter_lower_leg)
        struct_layout.addLayout(region_box)

        # Filter toggles — Type
        type_box = qt.QHBoxLayout()
        type_box.addWidget(qt.QLabel("Type:"))
        self.filter_bones = qt.QCheckBox("Bones")
        self.filter_muscles = qt.QCheckBox("Muscles")
        self.filter_skin = qt.QCheckBox("Skin")
        self.filter_landmarks = qt.QCheckBox("Landmarks")
        self.filter_bones.setChecked(True)
        self.filter_muscles.setChecked(True)
        self.filter_skin.setChecked(True)
        self.filter_landmarks.setChecked(True)
        type_box.addWidget(self.filter_bones)
        type_box.addWidget(self.filter_muscles)
        type_box.addWidget(self.filter_skin)
        type_box.addWidget(self.filter_landmarks)
        struct_layout.addLayout(type_box)

        # Structure checklist
        self.struct_list = qt.QListWidget()
        self.struct_list.setMaximumHeight(250)
        struct_layout.addWidget(self.struct_list)

        # --- Collapsible: Animation ---
        anim_section = ctk.ctkCollapsibleButton()
        anim_section.text = "Animation"
        self.layout.addWidget(anim_section)
        anim_layout = qt.QFormLayout(anim_section)

        # Play / Stop button
        self.play_button = qt.QPushButton("Play")
        self.play_button.enabled = False
        self.play_button.checkable = True
        anim_layout.addRow(self.play_button)

        # Speed slider
        self.speed_slider = ctk.ctkSliderWidget()
        self.speed_slider.singleStep = 0.01
        self.speed_slider.minimum = 0.01
        self.speed_slider.maximum = 0.50
        self.speed_slider.value = 0.05
        self.speed_slider.decimals = 2
        anim_layout.addRow("Step size:", self.speed_slider)

        # Interval (ms)
        self.interval_spin = qt.QSpinBox()
        self.interval_spin.minimum = 10
        self.interval_spin.maximum = 500
        self.interval_spin.value = 30
        self.interval_spin.suffix = " ms"
        anim_layout.addRow("Interval:", self.interval_spin)

        # Animation timer
        self._anim_timer = qt.QTimer()
        self._anim_direction = 1  # +1 or -1

        # Spacer
        self.layout.addStretch(1)

        # --- Connections ---
        self.load_button.connect("clicked(bool)", self.on_load_clicked)
        self.beta_slider.connect("valueChanged(double)", self.on_beta_changed)
        self.component_selector.connect(
            "currentIndexChanged(int)", self.on_component_changed
        )
        self.reset_button.connect("clicked(bool)", self.on_reset_clicked)
        self.range_min_spin.connect("valueChanged(double)", self.on_range_changed)
        self.range_max_spin.connect("valueChanged(double)", self.on_range_changed)
        self.color_selector.connect("currentIndexChanged(int)", self.on_color_changed)
        self.dim_opacity_spin.connect("valueChanged(int)", self.on_structure_visibility_changed)
        self.struct_list.connect("itemChanged(QListWidgetItem*)", self.on_structure_visibility_changed)
        self.select_all_btn.connect("clicked(bool)", self.on_select_all)
        self.select_none_btn.connect("clicked(bool)", self.on_select_none)
        self.filter_left.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_right.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_center.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_hip_thigh.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_lower_leg.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_bones.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_muscles.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_skin.connect("stateChanged(int)", self.on_apply_filters)
        self.filter_landmarks.connect("stateChanged(int)", self.on_apply_filters)
        self.play_button.connect("toggled(bool)", self.on_play_toggled)
        self._anim_timer.connect("timeout()", self.on_anim_tick)

    # ---- Initial path resolution ----

    @staticmethod
    def _resolve_initial_npz_path():
        """
        Resolve the initial NPZ path from (priority order):
          1. Slicer command-line:  --LLSSMExplorer-npz <path>
          2. Environment variable: LLSSM_EXPLORER_NPZ
          3. Slicer persistent setting: LLSSMExplorer/NpzPath (saved from last session)
        Returns path string or None.
        """
        import sys

        # 1. Command-line argument
        argv = sys.argv
        for i, arg in enumerate(argv):
            if arg == "--LLSSMExplorer-npz" and i + 1 < len(argv):
                candidate = argv[i + 1]
                if os.path.isfile(candidate):
                    return candidate

        # 2. Environment variable
        env_path = os.environ.get("LLSSM_EXPLORER_NPZ")
        if env_path and os.path.isfile(env_path):
            return env_path

        # 3. Slicer persistent setting (from previous session)
        settings = qt.QSettings()
        saved = settings.value("LLSSMExplorer/NpzPath", "")
        if saved and os.path.isfile(saved):
            return saved

        # 4. Built-in default (convenience fallback)
        _DEFAULT_NPZ = "llssm.npz"
        if os.path.isfile(_DEFAULT_NPZ):
            return _DEFAULT_NPZ

        return None

    @staticmethod
    def _resolve_initial_color_table():
        """
        Resolve the initial color table from (priority order):
          1. Slicer command-line:  --LLSSMExplorer-color <name>
          2. Environment variable: LLSSM_EXPLORER_COLOR
        Returns color table name string or None (use default).
        """
        import sys

        argv = sys.argv
        for i, arg in enumerate(argv):
            if arg == "--LLSSMExplorer-color" and i + 1 < len(argv):
                return argv[i + 1]

        env_color = os.environ.get("LLSSM_EXPLORER_COLOR")
        if env_color:
            return env_color

        return None

    @staticmethod
    def _resolve_initial_view():
        """
        Resolve the initial camera view from (priority order):
          1. Slicer command-line:  --LLSSMExplorer-view <direction>
          2. Environment variable: LLSSM_EXPLORER_VIEW
        Valid values: Anterior, Posterior, Left, Right, Superior, Inferior
        Returns string or None (use default).
        """
        import sys

        argv = sys.argv
        for i, arg in enumerate(argv):
            if arg == "--LLSSMExplorer-view" and i + 1 < len(argv):
                return argv[i + 1]

        env_view = os.environ.get("LLSSM_EXPLORER_VIEW")
        if env_view:
            return env_view

        return None

    # ---- Camera ----

    _VIEW_PRESETS = {
        "Anterior":  ctk.ctkAxesWidget.Anterior,
        "Posterior": ctk.ctkAxesWidget.Posterior,
        "Left":      ctk.ctkAxesWidget.Left,
        "Right":     ctk.ctkAxesWidget.Right,
        "Superior":  ctk.ctkAxesWidget.Superior,
        "Inferior":  ctk.ctkAxesWidget.Inferior,
    }

    @staticmethod
    def _apply_camera_view(threeDView, view_spec):
        """
        Apply camera orientation.

        view_spec can be:
          - Named preset: "Anterior", "Posterior", "Left", ...
          - Euler angles in degrees: "rx,ry,rz"  (e.g. "0,180,0")
            rx = elevation (pitch), ry = azimuth (yaw), rz = roll
        """
        import math

        # Try named preset first
        preset = LLSSMExplorerWidget._VIEW_PRESETS.get(view_spec)
        if preset is not None:
            threeDView.lookFromAxis(preset)
            return

        # Try x,y,z degree parsing
        parts = view_spec.replace(" ", "").split(",")
        if len(parts) == 3:
            try:
                elev_deg, azim_deg, roll_deg = (float(p) for p in parts)
            except ValueError:
                threeDView.lookFromAxis(ctk.ctkAxesWidget.Posterior)
                return

            renderer = threeDView.renderWindow().GetRenderers().GetFirstRenderer()
            if not renderer:
                return
            camera = renderer.GetActiveCamera()
            renderer.ResetCamera()

            # Compute camera position from spherical angles
            dist = camera.GetDistance()
            focal = list(camera.GetFocalPoint())
            elev = math.radians(elev_deg)
            azim = math.radians(azim_deg)

            # Spherical to Cartesian (RAS: X=Right, Y=Anterior, Z=Superior)
            x = dist * math.cos(elev) * math.sin(azim)
            y = dist * math.cos(elev) * math.cos(azim)
            z = dist * math.sin(elev)

            camera.SetPosition(
                focal[0] + x,
                focal[1] + y,
                focal[2] + z,
            )
            camera.SetViewUp(0, 0, 1)

            # Apply roll
            if roll_deg != 0:
                camera.Roll(roll_deg)

            renderer.ResetCameraClippingRange()
            return

        # Fallback
        threeDView.lookFromAxis(ctk.ctkAxesWidget.Posterior)

    # ---- Slots ----

    def on_load_clicked(self):
        npz_path = self.path_edit.currentPath
        if not npz_path or not os.path.isfile(npz_path):
            self.status_label.text = "ERROR: file not found"
            return

        try:
            self.logic.load_ssm(npz_path)
            color_name = self._color_table_names[self.color_selector.currentIndex]
            self.logic.create_model_node(color_table_name=color_name)
        except Exception as e:
            self.status_label.text = f"ERROR: {e}"
            import traceback
            traceback.print_exc()
            return

        # Populate component selector
        self.component_selector.clear()
        for i in range(self.logic.num_components):
            self.component_selector.addItem(f"PC{i+1}  (component {i})")
        self.component_selector.setCurrentIndex(0)

        # Enable controls
        self.component_selector.enabled = True
        self.beta_slider.enabled = True
        self.beta_slider.value = 0.0
        self.reset_button.enabled = True
        self.play_button.enabled = True

        V = self.logic.v_template.shape[0]
        K = self.logic.num_components
        self.status_label.text = f"Loaded: {V:,} vertices, {K} components"
        self._update_betas_display()

        # Populate structure checklist
        self._populate_structure_list()

        # Center the 3D view, set camera view, fix clipping
        threeDView = slicer.app.layoutManager().threeDWidget(0).threeDView()
        threeDView.resetFocalPoint()
        view_spec = self._resolve_initial_view() or "Posterior"
        self._apply_camera_view(threeDView, view_spec)
        # Extend the near clipping plane so zooming in doesn't clip the mesh
        viewNode = threeDView.mrmlViewNode()
        if viewNode:
            viewNode.SetNearPlaneScale(0.01)

    def _update_betas_display(self):
        """Update the active betas summary text."""
        if self.logic.current_beta is None:
            return
        parts = []
        for i, b in enumerate(self.logic.current_beta):
            if abs(b) > 1e-6:
                parts.append(f"PC{i+1}:{b:+.2f}")
        if parts:
            self.betas_display.plainText = "  ".join(parts)
        else:
            self.betas_display.plainText = "(mean shape)"

    def on_beta_changed(self, value):
        if self.logic.v_template is None:
            return
        comp = self.component_selector.currentIndex
        self.logic.set_beta(comp, value)
        self._update_betas_display()

    def on_range_changed(self, _value=None):
        self.beta_slider.minimum = self.range_min_spin.value
        self.beta_slider.maximum = self.range_max_spin.value

    def on_component_changed(self, index):
        if self.logic.v_template is None or index < 0:
            return
        self.beta_slider.value = self.logic.current_beta[index]

    def on_reset_clicked(self):
        if self.logic.v_template is None:
            return
        self.logic.reset_betas()
        idx = self.component_selector.currentIndex
        self.beta_slider.value = self.logic.current_beta[max(idx, 0)]
        self.status_label.text = "Reset to mean shape."
        self._update_betas_display()

    # ---- Display ----

    def on_color_changed(self, index):
        if self.logic.model_node is None:
            return
        name = self._color_table_names[index]
        try:
            color_node = slicer.util.getNode(name)
        except slicer.util.MRMLNodeNotFoundException:
            self.status_label.text = f"Color table '{name}' not found"
            return
        display_node = self.logic.model_node.GetDisplayNode()
        if display_node:
            display_node.SetAndObserveColorNodeID(color_node.GetID())

    # ---- Structures ----

    # Classification rules based on label names from the NPZ meta
    _BONE_KEYWORDS = {
        "pelvis", "femur", "tibia", "fibula", "patella",
        "sacrum", "foot_bone", "talus",
    }
    _LANDMARK_MIN_ID = 60

    # (Region classification is now geometry-based using vertex Z coordinates)

    def _classify_all_structures(self):
        """
        Classify every label using vertex geometry (data-driven).

        Side:   centroid X > threshold  → "right" (screen right)
                centroid X < -threshold → "left"  (screen left)
                otherwise              → "center"

        Region: centroid Z relative to knee_z midpoint
                above → "hip_thigh", below → "lower_leg"
                skin / spread structures → "other"

        Type:   label name keyword matching (bone vs muscle vs skin vs landmark)
        """
        v = self.logic.v_template        # (V, 3)
        labels = self.logic.vertex_labels # (V,)
        lid_to_name = self.logic.label_id_to_name

        # Pre-compute per-label centroids
        centroids = {}  # lid -> (cx, cy, cz)
        for lid in sorted(lid_to_name.keys()):
            mask = labels == lid
            if not np.any(mask):
                continue
            centroids[lid] = v[mask].mean(axis=0)  # (3,)

        # Determine side threshold: use the X spread of all centroids
        all_cx = np.array([c[0] for c in centroids.values()])
        x_median = float(np.median(all_cx))
        # Structures near the midline are "center"
        x_spread = float(np.std(all_cx))
        side_threshold = x_spread * 0.15  # within 15% of std → center

        # Determine region boundary: estimate knee Z as the boundary
        # between hip/thigh and lower-leg.
        # Use tibia/femur centroids if available, otherwise use Z median.
        femur_zs = []
        tibia_zs = []
        for lid, name in lid_to_name.items():
            if lid not in centroids:
                continue
            if "femur" in name:
                femur_zs.append(centroids[lid][2])
            elif "tibia" in name:
                tibia_zs.append(centroids[lid][2])
        if femur_zs and tibia_zs:
            knee_z = (np.mean(femur_zs) + np.mean(tibia_zs)) / 2.0
        else:
            all_cz = np.array([c[2] for c in centroids.values()])
            knee_z = float(np.median(all_cz))

        # Classify each structure
        tags_dict = {}
        for lid in sorted(lid_to_name.keys()):
            name = lid_to_name[lid]
            tags = {}

            # --- Side (geometry-based) ---
            if lid in centroids:
                cx = centroids[lid][0] - x_median
                if cx > side_threshold:
                    tags["side"] = "right"
                elif cx < -side_threshold:
                    tags["side"] = "left"
                else:
                    tags["side"] = "center"
            else:
                tags["side"] = "center"

            # --- Type (name-based — this is genuinely semantic) ---
            if lid >= self._LANDMARK_MIN_ID:
                tags["type"] = "landmark"
            elif name == "skin":
                tags["type"] = "skin"
            elif any(kw in name for kw in self._BONE_KEYWORDS):
                tags["type"] = "bone"
            else:
                tags["type"] = "muscle"

            # --- Region (geometry-based) ---
            if name == "skin":
                tags["region"] = "other"
            elif lid in centroids:
                cz = centroids[lid][2]
                if cz >= knee_z:
                    tags["region"] = "hip_thigh"
                else:
                    tags["region"] = "lower_leg"
            else:
                tags["region"] = "other"

            tags_dict[lid] = tags

        return tags_dict

    def _populate_structure_list(self):
        """Fill the structure checklist from loaded label_id_to_name."""
        self.struct_list.blockSignals(True)
        self.struct_list.clear()

        if not self.logic.label_id_to_name:
            self.struct_list.blockSignals(False)
            return

        self._struct_tags = self._classify_all_structures()

        for lid in sorted(self.logic.label_id_to_name.keys()):
            name = self.logic.label_id_to_name[lid]
            tags = self._struct_tags.get(lid, {})

            side_char = {"left": "L", "right": "R", "center": "C"}.get(
                tags.get("side", "center"), "?"
            )
            region_char = {"hip_thigh": "H", "lower_leg": "Lo", "other": "-"}.get(
                tags.get("region", "other"), "?"
            )
            item = qt.QListWidgetItem(
                f"[{lid}] {name}  ({side_char} {region_char})"
            )
            item.setFlags(item.flags() | qt.Qt.ItemIsUserCheckable)
            item.setCheckState(qt.Qt.Checked)
            item.setData(qt.Qt.UserRole, lid)
            self.struct_list.addItem(item)

        self.struct_list.blockSignals(False)

    def _get_checked_label_ids(self):
        """Return set of label IDs that are currently checked."""
        checked = set()
        for i in range(self.struct_list.count):
            item = self.struct_list.item(i)
            if item.checkState() == qt.Qt.Checked:
                checked.add(item.data(qt.Qt.UserRole))
        return checked

    def on_structure_visibility_changed(self, *args):
        if self.logic.model_node is None:
            return
        checked = self._get_checked_label_ids()
        dim_alpha = self.dim_opacity_spin.value
        color_name = self._color_table_names[self.color_selector.currentIndex]
        self.logic.update_structure_colors(checked, dim_alpha, color_name)

    def _set_all_check_state(self, state):
        self.struct_list.blockSignals(True)
        for i in range(self.struct_list.count):
            self.struct_list.item(i).setCheckState(state)
        self.struct_list.blockSignals(False)
        self.on_structure_visibility_changed()

    def on_select_all(self):
        self.struct_list.blockSignals(True)
        for i in range(self.struct_list.count):
            self.struct_list.item(i).setCheckState(qt.Qt.Checked)
        self.struct_list.blockSignals(False)
        # Restore original label-based coloring (no RGBA override)
        color_name = self._color_table_names[self.color_selector.currentIndex]
        self.logic.restore_label_colors(color_name)

    def on_select_none(self):
        self._set_all_check_state(qt.Qt.Unchecked)

    def on_apply_filters(self, *args):
        """
        Apply combinable filter toggles to the structure checklist.
        A structure is checked if it passes ALL active filters (AND logic):
          - Side: Left ON → allow left+center, Right ON → allow right+center
          - Region: Hip/Thigh ON → allow hip_thigh, Lower Leg ON → allow lower_leg
          - Type: Bones/Muscles/Skin/Landmarks toggles
        """
        if self.logic.model_node is None:
            return

        allow_left = self.filter_left.isChecked()
        allow_right = self.filter_right.isChecked()
        allow_center = self.filter_center.isChecked()
        allow_hip_thigh = self.filter_hip_thigh.isChecked()
        allow_lower_leg = self.filter_lower_leg.isChecked()
        allow_bones = self.filter_bones.isChecked()
        allow_muscles = self.filter_muscles.isChecked()
        allow_skin = self.filter_skin.isChecked()
        allow_landmarks = self.filter_landmarks.isChecked()

        # If ALL checkboxes in a category are OFF, skip that filter (show all)
        side_any = allow_left or allow_right or allow_center
        region_any = allow_hip_thigh or allow_lower_leg
        type_any = allow_bones or allow_muscles or allow_skin or allow_landmarks

        self.struct_list.blockSignals(True)

        for i in range(self.struct_list.count):
            item = self.struct_list.item(i)
            lid = item.data(qt.Qt.UserRole)
            tags = self._struct_tags.get(lid, {})

            # Side filter (skip if no side checkbox is checked)
            side = tags.get("side", "center")
            side_ok = (not side_any) or (
                (side == "center" and allow_center) or
                (side == "left" and allow_left) or
                (side == "right" and allow_right)
            )

            # Region filter (skip if no region checkbox is checked)
            # "other" (skin etc.) always passes — filtered by Type instead
            region = tags.get("region", "other")
            region_ok = (not region_any) or (
                (region == "other") or
                (region == "hip_thigh" and allow_hip_thigh) or
                (region == "lower_leg" and allow_lower_leg)
            )

            # Type filter (skip if no type checkbox is checked)
            stype = tags.get("type", "muscle")
            type_ok = (not type_any) or (
                (stype == "bone" and allow_bones) or
                (stype == "muscle" and allow_muscles) or
                (stype == "skin" and allow_skin) or
                (stype == "landmark" and allow_landmarks)
            )

            if side_ok and region_ok and type_ok:
                item.setCheckState(qt.Qt.Checked)
            else:
                item.setCheckState(qt.Qt.Unchecked)

        self.struct_list.blockSignals(False)
        self.on_structure_visibility_changed()

    # ---- Animation ----

    def on_play_toggled(self, checked):
        if checked:
            self.play_button.text = "Stop"
            self._anim_direction = 1
            self._anim_timer.start(self.interval_spin.value)
        else:
            self.play_button.text = "Play"
            self._anim_timer.stop()

    def on_anim_tick(self):
        step = self.speed_slider.value * self._anim_direction
        new_val = self.beta_slider.value + step

        # Bounce at boundaries
        if new_val >= self.beta_slider.maximum:
            new_val = self.beta_slider.maximum
            self._anim_direction = -1
        elif new_val <= self.beta_slider.minimum:
            new_val = self.beta_slider.minimum
            self._anim_direction = 1

        self.beta_slider.value = new_val
        # Update timer interval in case user changed it during playback
        self._anim_timer.setInterval(self.interval_spin.value)

    def cleanup(self):
        self._anim_timer.stop()


# ============================================================
# Logic (data + computation)
# ============================================================
class LLSSMExplorerLogic(ScriptedLoadableModuleLogic):

    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.v_template = None   # (V, 3) float32
        self.shapedirs = None    # (V, 3, K) float32
        self.faces_np = None     # (F, 3) int
        self.vertex_labels = None
        self.num_components = 0
        self.current_beta = None # (K,) float64

        self.polydata = None
        self.model_node = None

    # ---- Load ----

    def load_ssm(self, npz_path):
        data = np.load(npz_path, allow_pickle=True)
        self.v_template = data["v_template"].astype(np.float32)
        self.shapedirs = data["shapedirs"].astype(np.float32)
        self.faces_np = data["faces"].astype(np.int64)
        self.num_components = self.shapedirs.shape[2]
        self.current_beta = np.zeros(self.num_components, dtype=np.float64)

        if "point_data__vertex_labels" in data.files:
            self.vertex_labels = data["point_data__vertex_labels"]
        else:
            self.vertex_labels = None

        # Parse label_id_to_name from meta (for future use)
        self.label_id_to_name = {}
        if "__meta__" in data.files:
            import json
            try:
                meta = json.loads(str(data["__meta__"]))
                raw = meta.get("label_id_to_name", {})
                self.label_id_to_name = {
                    int(k): str(v) for k, v in raw.items()
                }
            except Exception:
                pass

    # ---- Deformation ----

    def compute_deformed(self):
        """v_template + einsum('k,vck->vc', beta, shapedirs)"""
        offsets = np.einsum("k,vck->vc", self.current_beta, self.shapedirs)
        return (self.v_template + offsets).astype(np.float32)

    def set_beta(self, component_index, value):
        self.current_beta[component_index] = value
        deformed = self.compute_deformed()
        self._update_points(deformed)

    def reset_betas(self):
        self.current_beta[:] = 0.0
        self._update_points(self.v_template.copy())

    # ---- VTK mesh ----

    def create_model_node(self, color_table_name="GenericAnatomyColors"):
        """Build vtkPolyData + MRMLModelNode from current LLSSM data."""

        # Remove previous node if exists
        if self.model_node is not None:
            slicer.mrmlScene.RemoveNode(self.model_node)
            self.model_node = None

        polydata = vtk.vtkPolyData()

        # Points
        pts = vtk.vtkPoints()
        vtk_arr = numpy_to_vtk(
            self.v_template.ravel(), deep=True, array_type=vtk.VTK_FLOAT
        )
        vtk_arr.SetNumberOfComponents(3)
        pts.SetData(vtk_arr)
        polydata.SetPoints(pts)

        # Faces
        num_faces = self.faces_np.shape[0]
        # Build VTK cell array: [3, v0, v1, v2, 3, v0, v1, v2, ...]
        cells_np = np.empty((num_faces, 4), dtype=np.int64)
        cells_np[:, 0] = 3
        cells_np[:, 1:] = self.faces_np
        vtk_cells = vtk.vtkCellArray()
        vtk_id_arr = numpy_to_vtk(
            cells_np.ravel(), deep=True, array_type=vtk.VTK_ID_TYPE
        )
        vtk_cells.SetCells(num_faces, vtk_id_arr)
        polydata.SetPolys(vtk_cells)

        # Vertex labels as scalars (colour by structure)
        if self.vertex_labels is not None:
            labels_vtk = numpy_to_vtk(
                self.vertex_labels.astype(np.float32), deep=True,
                array_type=vtk.VTK_FLOAT,
            )
            labels_vtk.SetName("Labels")
            polydata.GetPointData().AddArray(labels_vtk)
            polydata.GetPointData().SetActiveScalars("Labels")

        # Normals (for good shading)
        normals_filter = vtk.vtkPolyDataNormals()
        normals_filter.SetInputData(polydata)
        normals_filter.ComputePointNormalsOn()
        normals_filter.SplittingOff()
        normals_filter.Update()
        polydata = normals_filter.GetOutput()

        self.polydata = polydata

        # Create MRML model node
        model_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLModelNode", "LLSSM_Model"
        )
        model_node.SetAndObservePolyData(polydata)

        # Display node
        display_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLModelDisplayNode"
        )
        model_node.SetAndObserveDisplayNodeID(display_node.GetID())

        if self.vertex_labels is not None:
            display_node.SetScalarVisibility(True)
            display_node.SetActiveScalarName("Labels")
            try:
                color_node = slicer.util.getNode(color_table_name)
            except slicer.util.MRMLNodeNotFoundException:
                color_node = slicer.util.getNode("GenericAnatomyColors")
            display_node.SetAndObserveColorNodeID(color_node.GetID())
            display_node.SetScalarRangeFlag(
                slicer.vtkMRMLDisplayNode.UseDataScalarRange
            )
        else:
            display_node.SetColor(0.8, 0.8, 0.9)

        display_node.SetVisibility(True)

        # Widen camera clipping range so the model doesn't disappear on zoom
        threeDWidget = slicer.app.layoutManager().threeDWidget(0)
        if threeDWidget:
            renderer = threeDWidget.threeDView().renderWindow().GetRenderers().GetFirstRenderer()
            if renderer:
                camera = renderer.GetActiveCamera()
                camera.SetClippingRange(0.1, 100000)
                renderer.ResetCameraClippingRange()

        self.model_node = model_node

    def update_structure_colors(self, selected_label_ids, dim_alpha=40,
                                color_table_name="GenericAnatomyColors"):
        """
        Set per-vertex RGBA: selected structures get full color,
        unselected get dim grey with low alpha.
        """
        if self.polydata is None or self.vertex_labels is None:
            return

        # Get color lookup table from Slicer
        try:
            color_node = slicer.util.getNode(color_table_name)
        except slicer.util.MRMLNodeNotFoundException:
            color_node = slicer.util.getNode("GenericAnatomyColors")
        lut = color_node.GetLookupTable()

        V = len(self.vertex_labels)
        rgba = np.zeros((V, 4), dtype=np.uint8)

        # Dim color for unselected
        rgba[:, 0] = 180  # R
        rgba[:, 1] = 180  # G
        rgba[:, 2] = 180  # B
        rgba[:, 3] = dim_alpha

        # Vectorized: get colors for all unique labels
        num_lut = lut.GetNumberOfTableValues()
        unique_labels = np.unique(self.vertex_labels)
        color4 = [0.0, 0.0, 0.0, 1.0]
        for lid in unique_labels:
            lid_int = int(lid)
            mask = self.vertex_labels == lid
            if lid_int in selected_label_ids:
                # Label 0 is often "background=black" in most color tables;
                # use a skin-like color instead
                if lid_int == 0:
                    rgba[mask, 0] = 255
                    rgba[mask, 1] = 224
                    rgba[mask, 2] = 189
                    rgba[mask, 3] = 255
                else:
                    lut.GetTableValue(lid_int % num_lut, color4)
                    rgba[mask, 0] = int(color4[0] * 255)
                    rgba[mask, 1] = int(color4[1] * 255)
                    rgba[mask, 2] = int(color4[2] * 255)
                    rgba[mask, 3] = 255

        # Apply to polydata
        vtk_colors = numpy_to_vtk(
            rgba.ravel(), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR
        )
        vtk_colors.SetNumberOfComponents(4)
        vtk_colors.SetName("RGBA")
        self.polydata.GetPointData().SetScalars(vtk_colors)
        self.polydata.Modified()

        # Switch display to use direct colors (not lookup table)
        display_node = self.model_node.GetDisplayNode()
        if display_node:
            display_node.SetScalarVisibility(False)
            display_node.SetScalarVisibility(True)
            display_node.SetActiveScalarName("RGBA")
            display_node.SetScalarRangeFlag(
                slicer.vtkMRMLDisplayNode.UseDirectMapping
            )

        if self.model_node is not None:
            self.model_node.Modified()

    def restore_label_colors(self, color_table_name="GenericAnatomyColors"):
        """Switch back to label-based coloring (undo structure highlight)."""
        if self.polydata is None or self.vertex_labels is None:
            return

        self.polydata.GetPointData().SetActiveScalars("Labels")
        self.polydata.Modified()

        display_node = self.model_node.GetDisplayNode()
        if display_node:
            try:
                color_node = slicer.util.getNode(color_table_name)
            except slicer.util.MRMLNodeNotFoundException:
                color_node = slicer.util.getNode("GenericAnatomyColors")
            display_node.SetAndObserveColorNodeID(color_node.GetID())
            display_node.SetActiveScalarName("Labels")
            display_node.SetScalarRangeFlag(
                slicer.vtkMRMLDisplayNode.UseDataScalarRange
            )

        if self.model_node is not None:
            self.model_node.Modified()

    def _update_points(self, vertices):
        """Fast in-place update of vtkPolyData points."""
        if self.polydata is None:
            return

        vtk_arr = numpy_to_vtk(
            vertices.ravel(), deep=True, array_type=vtk.VTK_FLOAT
        )
        vtk_arr.SetNumberOfComponents(3)
        self.polydata.GetPoints().SetData(vtk_arr)
        self.polydata.GetPoints().Modified()
        self.polydata.Modified()

        if self.model_node is not None:
            self.model_node.Modified()
