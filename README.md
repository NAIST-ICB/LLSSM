# Auto-FWLL-SSM (LLSSM)

**Automated Female Whole Lower Limb Statistical Shape Model**

Official implementation of the paper:

> **Multi-structure statistical shape model of the lower limb musculoskeletal system in patients with hip osteoarthritis**
>
> Yuto Masaki, Yoshito Otake, Mazen Soufi, Keisuke Uemura, Masaki Takao, Takuma Miyamoto, Yasuhito Tanaka, Seiji Okada, Nobuhiko Sugano, Yoshinobu Sato
>
> *Under review*

## Overview

Auto-FWLL-SSM is a Statistical Shape Model (SSM) of the entire lower limb musculoskeletal system, constructed from CT scans of 560 Japanese female patients with hip osteoarthritis (HOA). The model integrates 60 anatomical structures including bones, muscles, and skin using deep learning-based automated segmentation and PCA-based shape modeling.

This repository provides tools for:

- Reading and writing multi-structure polygon meshes (VTK format)
- Decomposing labeled meshes into individual anatomical structures
- Applying SSM-based shape deformation using PCA components

## Requirements

- Python 3.x
- NumPy
- VTK

### Optional

- PyTorch (for GPU/batched SSM operations)
- SimpleITK (for reading `.mha`/`.mhd`/`.nii.gz` files)

## SSM Models

> Under construction

## Usage

<!--
```bash
# Read a labeled mesh
python multi_structure_mesh_reader.py <input.vtk> <output.vtk> <point_data_name>
```
-->


```bash
# Read a labeled mesh and export with structure decomposition
python multi_structure_mesh_reader.py <input.vtk> <output.vtk> <point_data_name>

# Apply SSM deformation
python multi_structure_mesh_reader.py <input.vtk> <output.vtk> <point_data_name> \
    <mean_mesh.vtk> <components.txt> <variances.txt> <num_components>

# Apply SSM deformation with explicit shape parameters
python multi_structure_mesh_reader.py <input.vtk> <output.vtk> <point_data_name> \
    <mean_mesh.vtk> <components.txt> <variances.txt> <num_components> <shape_params.txt>
```


## Simple Visualization

The mean shape of the Auto-FWLL-SSM, visualized in 3D Slicer with per-structure color labels:

<!-- 
![Mean shape visualized in 3D Slicer](fig/mean_shape_3dslicer.png) 
-->

`visualize_3dslicer.py` is a 3D Slicer Python script for automated mesh visualization and screen capture.

```bash
# Usage (run via 3D Slicer's built-in Python)
"<path/to/Slicer>" --python-script visualize_3dslicer.py <input.vtk> \
    --capture <output.png> --auto-crop
```

Options: `--array`, `--colormap`, `--opacity`, `--range-min`, `--range-max`, `--rotx/y/z`, `--zoom`, `--pan-x/y`, `--parallel`, `--fov`, `--capture`.

## LLSSM Explorer — 3D Slicer Plugin

A 3D Slicer plugin (`LLSSMExplorer/`) for real-time interactive exploration of PCA shape modes in the LLSSM. Load an `.npz` model file and interactively deform the mean shape along each principal component using a slider (-3σ to +3σ).

Features:
- Real-time PCA shape deformation with cumulative multi-component support
- Anatomical structure filtering by side (Left/Right/Center), region (Hip-Thigh/Lower Leg), and type (Bones/Muscles/Skin/Landmarks)
- Animation playback for automatic mode oscillation
- Configurable color tables and structure highlighting

```bat
"<path/to/Slicer>" ^
    --additional-module-paths "path\to\LLSSMExplorer" ^
    --LLSSMExplorer-npz "path\to\model.npz" ^
    --LLSSMExplorer-view "Posterior"
```

<!-- See [`LLSSMExplorer/README.md`](LLSSMExplorer/README.md) for full documentation. -->

## Data Availability

> Under construction

<!-- The Auto-FWLL-SSM data (mean mesh, PCA components, and variances) will be made publicly available upon publication. -->

## Citation

> Under construction

<!--
If you use this code or data in your research, please cite:

```bibtex
@article{masaki2026autofwllssm,
  title={Multi-structure statistical shape model of the lower limb musculoskeletal system in patients with hip osteoarthritis},
  author={Masaki, Yuto and Otake, Yoshito and Soufi, Mazen and Uemura, Keisuke and Takao, Masaki and Miyamoto, Takuma and Tanaka, Yasuhito and Okada, Seiji and Sugano, Nobuhiko and Sato, Yoshinobu},
  journal={},
  year={}
}
```
-->

