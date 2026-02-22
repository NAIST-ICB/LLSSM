import sys
import os
import numpy as np
from collections import OrderedDict

import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa

torch = None
try:
    import torch
except Exception:
    print("exception : import torch")

sitk = None
try:
    import SimpleITK as sitk
except Exception:
    print("exception : import SimpleITK as sitk")

from resource import _name2id
from resource import _strcture_name_list_without_skin
from resource import _strcture_name_list

_default_point_data_name = "label"

_NUMPY_DTYPE_TO_VTK_TYPE = {
    np.float64: "double",
    np.float32: "float",
    np.int32: "int",
    np.uint32: "unsigned_int",
    np.int8: "char",
    np.uint8: "unsigned_char",
    np.int16: "short",
    np.uint16: "unsigned_short",
    np.int64: "long",
    np.uint64: "unsigned_long",
}


def composite_path(filename):
    dirname, basefile = os.path.split(filename)
    basename, ext = os.path.splitext(basefile)
    return dirname, basename, ext


def check_file(in_file):
    return (isinstance(in_file, str)
            and os.path.exists(in_file)
            and os.path.isfile(in_file))


def check_folder(in_folder):
    return (isinstance(in_folder, str)
            and os.path.exists(in_folder)
            and os.path.isdir(in_folder))


def read_matrix(in_file, dtype=None):
    if dtype is None:
        dtype = "float"
    matrix = None
    if check_file(in_file):
        _, _, in_ext = composite_path(filename=in_file)
        if in_ext in (".txt", ".dat"):
            matrix = np.loadtxt(fname=in_file, dtype=dtype, comments="#",
                                delimiter=None, converters=None, skiprows=0,
                                usecols=None, unpack=False, ndmin=0)
        elif in_ext == ".npy":
            matrix = np.load(in_file)
    return matrix


def read_matrix_sitk(in_file, dtype=None):
    matrix = read_matrix(in_file=in_file, dtype=dtype)
    if matrix is None and check_file(in_file):
        _, _, in_ext = composite_path(filename=in_file)
        if in_ext in (".mha", ".mhd", ".gz", ".nii.gz"):
            image_sitk = sitk.ReadImage(in_file)
            matrix = sitk.GetArrayFromImage(image_sitk)
    return matrix


def write_matrix(out_file, x, fmt=None, do_calc_fmt_in=None):
    if not isinstance(out_file, str):
        return
    out_folder, out_fname, out_ext = composite_path(out_file)

    if isinstance(x, dict):
        out_npz_file = os.path.join(out_folder, out_fname + ".npz")
        np.savez(out_npz_file, **x)
    else:
        if do_calc_fmt_in:
            x_dtype = x.dtype
            if x_dtype in (np.float16, np.float32, np.float64):
                fmt = "%lf"
            else:
                fmt = "%d"
        if fmt is None:
            fmt = "%lf"
        if out_ext == ".txt":
            np.savetxt(fname=out_file, X=x, fmt=fmt, delimiter="\t",
                       newline="\n", header="", footer="", comments="#")
        elif out_ext == ".npy":
            np.save(out_file, x)


class Mesh():
    def __init__(self, points, faces, point_data=None, point_data_dict=None, point_data_name=None):
        self._points = points
        self._faces = faces
        self._point_data = point_data
        self._point_data_dict = point_data_dict
        self._point_data_name = point_data_name

    @property
    def points(self):
        return self._points

    @property
    def faces(self):
        return self._faces

    @property
    def point_data(self):
        return self._point_data

    @property
    def point_data_name(self):
        return self._point_data_name

    @property
    def point_data_dict(self):
        return self._point_data_dict


class MeshReader():
    def __init__(self, in_file, point_data_name=None, point_data_name_list=None,
                 should_read_point_data=None, should_read_point_data_list=None):
        self._in_file = in_file
        self._point_data_name = point_data_name
        self._point_data_name_list = point_data_name_list
        self._should_read_point_data = should_read_point_data
        self._should_read_point_data_list = should_read_point_data_list
        self._mesh = None

    @property
    def in_file(self):
        return self._in_file

    @property
    def point_data_name(self):
        return self._point_data_name

    @property
    def point_data_name_list(self):
        return self._point_data_name_list

    @property
    def should_read_point_data(self):
        return self._should_read_point_data

    @property
    def should_read_point_data_list(self):
        return self._should_read_point_data_list

    @property
    def mesh(self):
        return self._mesh

    def run(self):
        self._mesh = self.read(in_file=self.in_file,
                               point_data_name=self.point_data_name,
                               point_data_name_list=self.point_data_name_list)
        return self.mesh

    def read(self, in_file, point_data_name, point_data_name_list=None):
        vtk_polydata = self.read_polygon(in_file=in_file)
        if vtk_polydata is None:
            return None
        points = self.get_points(vtk_polydata=vtk_polydata)
        faces = self.get_faces(vtk_polydata=vtk_polydata)

        vtk_polydata_dsa = None
        if isinstance(point_data_name, str) or isinstance(point_data_name_list, (list, tuple, dict, str)):
            vtk_polydata_dsa = dsa.WrapDataObject(vtk_polydata)

        point_data = None
        point_data_dict = None
        if vtk_polydata_dsa is not None:
            if self.should_read_point_data and isinstance(point_data_name, str):
                point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=point_data_name)

            if self.should_read_point_data_list and point_data_name_list is not None:
                if isinstance(point_data_name_list, str):
                    point_data_raw = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=point_data_name_list)
                    if point_data_raw is not None:
                        point_data_dict = OrderedDict()
                        point_data_dict[point_data_name_list] = point_data_raw
                elif isinstance(point_data_name_list, (list, tuple)):
                    point_data_dict = self.read_point_data_list(vtk_polydata_dsa=vtk_polydata_dsa, name_list=point_data_name_list)
                elif isinstance(point_data_name_list, dict):
                    point_data_dict = self.read_point_data_dict(vtk_polydata_dsa=vtk_polydata_dsa, name_dict=point_data_name_list)

        return Mesh(points=points, faces=faces, point_data=point_data,
                    point_data_dict=point_data_dict, point_data_name=point_data_name)

    def read_polygon(self, in_file):
        _, _, in_ext = composite_path(filename=in_file)

        _READERS = {
            '.stl': vtk.vtkSTLReader,
            '.ply': vtk.vtkPLYReader,
            '.obj': vtk.vtkOBJReader,
            '.vtp': vtk.vtkXMLPolyDataReader,
            '.vtk': vtk.vtkPolyDataReader,
        }

        reader_cls = _READERS.get(in_ext)
        if reader_cls is None:
            return None

        reader = reader_cls()
        reader.SetFileName(in_file)
        reader.Update()
        return reader.GetOutput()

    def get_points(self, vtk_polydata):
        vtk_points = vtk_polydata.GetPoints()
        vtkarray_points = vtk_points.GetData()
        return vtk_to_numpy(vtkarray_points)

    def get_faces(self, vtk_polydata):
        vtk_cells = vtk_polydata.GetPolys()
        n_cells = vtk_cells.GetNumberOfCells()
        vtkarray_cells = vtk_cells.GetData()
        n_cols = vtkarray_cells.GetNumberOfValues() // n_cells
        cells = vtk_to_numpy(vtkarray_cells)
        cells = cells.reshape((-1, n_cols))
        return cells[:, 1:]

    def read_point_data(self, vtk_polydata_dsa, name):
        if check_file(name):
            return read_matrix_sitk(in_file=name, dtype=None)
        return vtk_polydata_dsa.PointData[name]

    def read_point_data_list(self, vtk_polydata_dsa, name_list):
        point_data_dict = None
        for i in range(len(name_list)):
            name = name_list[i]
            if name is None:
                continue
            point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=name)
            if point_data is None:
                continue
            if point_data_dict is None:
                point_data_dict = OrderedDict()
            point_data_dict[name] = point_data
        return point_data_dict

    def read_point_data_dict(self, vtk_polydata_dsa, name_dict):
        point_data_dict = None
        for k in name_dict.keys():
            if k is None:
                continue
            v = name_dict.get(k)
            if v is None:
                v = k
            point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=v)
            if point_data is None:
                point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=k)
            if point_data is None:
                continue
            if point_data_dict is None:
                point_data_dict = OrderedDict()
            point_data_dict[k] = point_data
        return point_data_dict


class MeshDisassembler():
    def __init__(self, mesh, structure_id_list, structure_name_list=None, name2id=None):
        self._mesh = mesh
        self._structure_id_list = structure_id_list
        self._structure_name_list = structure_name_list
        self._name2id = name2id
        self._disassembled_mesh = None

        if self._structure_id_list is None:
            if isinstance(structure_name_list, (list, tuple)) and isinstance(name2id, dict):
                self._structure_id_list = self.make_structure_id_list(
                    structure_name_list=structure_name_list, name2id=name2id)

    @property
    def mesh(self):
        return self._mesh

    @property
    def structure_id_list(self):
        return self._structure_id_list

    @property
    def disassembled_mesh(self):
        return self._disassembled_mesh

    def run(self):
        self._disassembled_mesh = self.process(mesh=self.mesh, structure_id_list=self.structure_id_list)
        return self.disassembled_mesh

    def process(self, mesh, structure_id_list):
        if structure_id_list is None:
            return None
        if not isinstance(structure_id_list, (np.ndarray, tuple, list, dict)):
            return None

        if isinstance(structure_id_list, dict):
            structure_id_list = sorted(set(structure_id_list.values()))

        return self.assemble(mesh=mesh, structure_id_list=structure_id_list)

    def disassemble(self, structure_id):
        self._disassembled_mesh = self.decompose(mesh=self.mesh, structure_id=structure_id)

    def assemble(self, mesh, structure_id_list):
        if mesh is None or structure_id_list is None:
            return None

        mesh_points_max = 0
        points_list = []
        faces_list = []
        point_data_list = []
        point_data_dict_list = None

        for structure_id in structure_id_list:
            decomposed_mesh = self.decompose(mesh=mesh, structure_id=structure_id)
            if decomposed_mesh is None:
                continue

            len_of_decomposed_mesh_points = len(decomposed_mesh.points)
            decomposed_mesh_faces_updated = decomposed_mesh.faces + mesh_points_max

            points_list.append(decomposed_mesh.points)
            faces_list.append(decomposed_mesh_faces_updated)
            point_data_list.append(decomposed_mesh.point_data)

            if isinstance(decomposed_mesh.point_data_dict, dict):
                for k in mesh.point_data_dict.keys():
                    if k is None:
                        continue
                    v = decomposed_mesh.point_data_dict.get(k)
                    if v is None:
                        continue
                    if point_data_dict_list is None:
                        point_data_dict_list = OrderedDict()
                    if point_data_dict_list.get(k) is None:
                        point_data_dict_list[k] = []
                    point_data_dict_list[k].append(v)

            mesh_points_max += len_of_decomposed_mesh_points

        points_concat = np.concatenate(points_list, axis=0) if points_list else None
        faces_concat = np.concatenate(faces_list, axis=0) if faces_list else None
        point_data_concat = np.concatenate(point_data_list, axis=0) if point_data_list else None

        point_data_dict_concat = None
        if isinstance(point_data_dict_list, dict):
            for k, v in point_data_dict_list.items():
                if k is None or v is None or len(v) == 0:
                    continue
                if point_data_dict_concat is None:
                    point_data_dict_concat = OrderedDict()
                point_data_dict_concat[k] = np.concatenate(v, axis=0)

        return Mesh(points=points_concat, faces=faces_concat,
                    point_data=point_data_concat, point_data_dict=point_data_dict_concat,
                    point_data_name=mesh.point_data_name)

    def decompose(self, mesh, structure_id):
        if mesh is None or mesh.point_data is None:
            return None

        area = (mesh.point_data == structure_id)
        if np.count_nonzero(area) <= 0:
            return None

        indices = np.where(area)
        pts = mesh.points[area]
        pd = mesh.point_data[area] if mesh.point_data is not None else None

        indices_dict = self.convert_list_to_dict(list_like_object=indices[0], do_reverse_key_value=True)
        fcs = self.extract_faces(valid_vertex_id_table=indices_dict, faces=mesh.faces)

        pd_dict = None
        if isinstance(mesh.point_data_dict, dict):
            for k, v in mesh.point_data_dict.items():
                if k is None or v is None:
                    continue
                vpd = v[area]
                if vpd is None:
                    continue
                if pd_dict is None:
                    pd_dict = OrderedDict()
                pd_dict[k] = vpd

        return Mesh(points=pts, faces=fcs, point_data=pd,
                    point_data_dict=pd_dict, point_data_name=mesh.point_data_name)

    def convert_list_to_dict(self, list_like_object, do_reverse_key_value=None):
        if do_reverse_key_value is None:
            do_reverse_key_value = True
        dict_object = OrderedDict()
        for i in range(len(list_like_object)):
            value = list_like_object[i]
            if do_reverse_key_value:
                dict_object[value] = i
            else:
                dict_object[i] = value
        return dict_object

    def convert_list_to_dict_tuple0(self, tuple_list_like_object, do_reverse_key_value=None):
        return self.convert_list_to_dict(list_like_object=tuple_list_like_object[0],
                                         do_reverse_key_value=do_reverse_key_value)

    def extract_faces(self, valid_vertex_id_table, faces):
        new_face_list = []
        for face in faces:
            face_0 = valid_vertex_id_table.get(face[0])
            face_1 = valid_vertex_id_table.get(face[1])
            face_2 = valid_vertex_id_table.get(face[2])
            if face_0 is not None and face_1 is not None and face_2 is not None:
                new_face_list.append(np.array((face_0, face_1, face_2)))
        return np.array(new_face_list)

    def make_structure_id_list(self, structure_name_list, name2id):
        structure_id_list = []
        for structure_name in structure_name_list:
            if structure_name is None:
                continue
            structure_id = name2id.get(structure_name)
            if structure_id is None:
                continue
            structure_id_list.append(structure_id)
        return structure_id_list if structure_id_list else None


class MeshVTKWriter():
    def __init__(self, mesh, out_file):
        self._mesh = mesh
        self._out_file = out_file

    @property
    def mesh(self):
        return self._mesh

    @property
    def out_file(self):
        return self._out_file

    def write_vtk_polygon_pts(self, ff, points):
        if points is None:
            return
        ff.write("POINTS {} float\n".format(len(points)))
        for point in points:
            ff.write("{} {} {} \n".format(point[0], point[1], point[2]))

    def write_vtk_polygon_faces(self, ff, faces, mesh_ndim=3):
        if faces is None:
            return
        faces_shape = faces.shape
        len_of_faces = faces_shape[0]
        dim_of_content = (faces_shape[0] * faces_shape[1]) + len_of_faces
        ff.write("POLYGONS {} {}\n".format(len_of_faces, dim_of_content))
        for face in faces:
            ff.write("{} {} {} {} \n".format(mesh_ndim, face[0], face[1], face[2]))

    def write_vtk_polygon_scalar_vector(self, ff, name, scalar_vector):
        if scalar_vector is None:
            return
        if name is None:
            name = "scalar"

        len_of_scalar_vector = len(scalar_vector)
        vtk_type = _NUMPY_DTYPE_TO_VTK_TYPE.get(scalar_vector.dtype.type, "double")
        ff.write("{} 1 {} {}\n".format(name, len_of_scalar_vector, vtk_type))

        for scalar in scalar_vector:
            try:
                len_of_scalar = len(scalar)
                text = " ".join(str(scalar[j]) for j in range(len_of_scalar)) + "\n"
                ff.write(text)
            except TypeError:
                ff.write("{} \n".format(scalar))

    def _write_vtk_polygon_file(self, filename, points, faces, scalar_vector=None,
                                scalar_vector_name=None, scalar_vector_dict=None,
                                do_cast_points=None, mesh_ndim=None, write_mode=None):
        if filename is None or points is None or faces is None:
            return False

        if do_cast_points is None:
            do_cast_points = True
        if mesh_ndim is None:
            mesh_ndim = points.shape[1]
        if write_mode is None:
            write_mode = "w"

        if do_cast_points:
            points = points.astype(np.float32)

        len_of_points = len(points)

        with open(filename, write_mode) as ff:
            ff.write("# vtk DataFile Version 4.2\n")
            ff.write("vtk output\n")
            ff.write("BINARY\n" if write_mode == "wb" else "ASCII\n")
            ff.write("DATASET POLYDATA\n")
            self.write_vtk_polygon_pts(ff=ff, points=points)
            self.write_vtk_polygon_faces(ff=ff, faces=faces, mesh_ndim=mesh_ndim)

            if scalar_vector is not None or isinstance(scalar_vector_dict, dict):
                ff.write("POINT_DATA {}\n".format(len_of_points))
                counter = 0
                if scalar_vector is not None:
                    counter += 1
                if isinstance(scalar_vector_dict, dict):
                    counter += len(scalar_vector_dict)
                ff.write("FIELD FieldData {}\n".format(counter))

                if scalar_vector is not None:
                    self.write_vtk_polygon_scalar_vector(
                        ff=ff, name=scalar_vector_name, scalar_vector=scalar_vector)

                if isinstance(scalar_vector_dict, dict):
                    for key, value in scalar_vector_dict.items():
                        if value is not None and len_of_points == len(value):
                            self.write_vtk_polygon_scalar_vector(
                                ff=ff, name=key, scalar_vector=value)

        return True

    def write_vtk_polygon_file(self, do_cast_points=None, mesh_ndim=None, write_mode=None):
        out_file = self.out_file
        out_dirname, _, _ = composite_path(filename=out_file)
        os.makedirs(out_dirname, exist_ok=True)
        return self._write_vtk_polygon_file(
            filename=out_file, points=self.mesh.points, faces=self.mesh.faces,
            scalar_vector=self.mesh.point_data, scalar_vector_name=self.mesh.point_data_name,
            scalar_vector_dict=self.mesh.point_data_dict,
            do_cast_points=do_cast_points, mesh_ndim=mesh_ndim, write_mode=write_mode)

    def run(self):
        return self.write_vtk_polygon_file()


class SSM():
    def __init__(self, mean_points, components, standard_deviations,
                 mean_points_flatten=None, number_of_batch=None, logger=None):
        self._mean_points = mean_points
        if mean_points_flatten is None and mean_points is not None:
            if isinstance(mean_points, np.ndarray):
                mean_points_flatten = mean_points.ravel()
            else:
                mean_points_flatten = torch.flatten(mean_points)
        self._mean_points_flatten = mean_points_flatten
        self._components = components
        self._standard_deviations = standard_deviations

        if number_of_batch is None:
            number_of_batch = 1
        self._number_of_batch = number_of_batch
        self._mean_points_batch = None
        self._mean_points_batch_flatten = None
        self._logger = logger

    @property
    def mean_points(self):
        return self._mean_points

    @property
    def components(self):
        return self._components

    @property
    def standard_deviations(self):
        return self._standard_deviations

    @property
    def number_of_batch(self):
        return self._number_of_batch

    @property
    def mean_points_flatten(self):
        return self._mean_points_flatten

    def _deform_term(self, shape_parameter, standard_deviations=None):
        if standard_deviations is None:
            standard_deviations = self.standard_deviations
        shape_score = standard_deviations * shape_parameter
        if isinstance(self.components, np.ndarray):
            return np.matmul(shape_score, self.components)
        else:
            return torch.matmul(shape_score, self.components)

    def _deform(self, shape_parameter, standard_deviations=None, mean_points_flatten=None):
        if mean_points_flatten is None:
            mean_points_flatten = self.mean_points_flatten
        deform_term = self._deform_term(shape_parameter=shape_parameter,
                                        standard_deviations=standard_deviations)
        return mean_points_flatten + deform_term

    def deform(self, shape_parameter):
        if shape_parameter is None:
            return self.mean_points
        deformed_verts_flatten = self._deform(
            shape_parameter=shape_parameter,
            standard_deviations=self.standard_deviations,
            mean_points_flatten=self.mean_points_flatten)
        return deformed_verts_flatten.reshape(-1, self.mean_points.shape[1])

    def deform_batch(self, shape_parameter):
        if shape_parameter is None:
            return self._mean_points_batch
        if isinstance(self.standard_deviations, np.ndarray):
            standard_deviations = np.expand_dims(self.standard_deviations, axis=0)
        else:
            standard_deviations = torch.unsqueeze(self.standard_deviations, dim=0)
        deformed_verts_flatten = self._deform(
            shape_parameter=shape_parameter,
            standard_deviations=standard_deviations,
            mean_points_flatten=self.mean_points_flatten)
        return deformed_verts_flatten.reshape(self._number_of_batch, -1, self.mean_points.shape[1])

    def calculate_shape_parameter(self, points):
        points_flatten = points.ravel()
        diff = points_flatten - self.mean_points_flatten
        score = np.matmul(self.components, diff)
        return score / self.standard_deviations


def generate_deformed_meshes(shape_parameter, ssm_mean_mesh, ssm_calculator,
                             structure_id_list, structure_name_list=None, name2id=None):
    deformed_points = ssm_calculator.deform(shape_parameter=shape_parameter)
    deformed_ssm_mesh = Mesh(
        points=deformed_points,
        faces=ssm_mean_mesh.faces,
        point_data=ssm_mean_mesh.point_data,
        point_data_dict=ssm_mean_mesh.point_data_dict,
        point_data_name=ssm_mean_mesh.point_data_name)

    deformed_ssm_mesh_disassembler = MeshDisassembler(
        mesh=deformed_ssm_mesh,
        structure_id_list=structure_id_list,
        structure_name_list=structure_name_list,
        name2id=name2id)
    disassembled_deformed_ssm_mesh = deformed_ssm_mesh_disassembler.run()
    return deformed_ssm_mesh, disassembled_deformed_ssm_mesh


def main(
    sys_argv,
    logger=None,
    mode=None,
    input_file=None,
    output_file=None,
    point_data_name=None,
    point_data_name_list=None,
    should_read_point_data=None,
    should_read_point_data_list=None,
    ssm_mean_mesh_file=None,
    ssm_components_file=None,
    ssm_variances_file=None,
    ssm_point_data_name=None,
    ssm_point_data_name_list=None,
    ssm_should_read_point_data=None,
    ssm_should_read_point_data_list=None,
    number_of_components=None,
    shape_parameter_file=None,
    structure_id_list=None,
    structure_name_list=None,
    name2id=None,
    do_overwrite=None
    ):

    if sys_argv is not None and isinstance(sys_argv, (list, tuple)):
        len_of_sys_argv = len(sys_argv)
        if len_of_sys_argv > 0:
            if mode is None:
                mode = 0
            if mode == 0:
                if len_of_sys_argv == 2:
                    input_file = str(sys_argv[1])
                elif len_of_sys_argv == 3:
                    input_file = str(sys_argv[1])
                    output_file = str(sys_argv[2])
                elif len_of_sys_argv == 4:
                    input_file = str(sys_argv[1])
                    output_file = str(sys_argv[2])
                    point_data_name = str(sys_argv[3])
                elif len_of_sys_argv == 8:
                    input_file = str(sys_argv[1])
                    output_file = str(sys_argv[2])
                    point_data_name = str(sys_argv[3])
                    ssm_mean_mesh_file = str(sys_argv[4])
                    ssm_components_file = str(sys_argv[5])
                    ssm_variances_file = str(sys_argv[6])
                    number_of_components = str(sys_argv[7])
                elif len_of_sys_argv == 9:
                    input_file = str(sys_argv[1])
                    output_file = str(sys_argv[2])
                    point_data_name = str(sys_argv[3])
                    ssm_mean_mesh_file = str(sys_argv[4])
                    ssm_components_file = str(sys_argv[5])
                    ssm_variances_file = str(sys_argv[6])
                    number_of_components = str(sys_argv[7])
                    shape_parameter_file = str(sys_argv[8])
                else:
                    print("1:input_file(str)[.vtk]")
                    print("2:output_file(str)[.vtk]")
                    print("3:point_data_name(str)")
                    print("4:ssm_mean_mesh_file(str)[.vtk]")
                    print("5:ssm_components_file(str)[.txt,.mhd/.mha]")
                    print("6:ssm_variances_file(str)[.txt]")
                    print("7:number_of_components(str,int)", "_none_")
                    print("8:shape_parameter_file(str)[.txt]")
                    return

    print("input_file", input_file)
    print("output_file", output_file)
    print("point_data_name", point_data_name)
    print("ssm_mean_mesh_file", ssm_mean_mesh_file)
    print("ssm_components_file", ssm_components_file)
    print("ssm_variances_file", ssm_variances_file)
    print("number_of_components", number_of_components)
    print("shape_parameter_file", shape_parameter_file)

    if point_data_name is None:
        point_data_name = _default_point_data_name
    if ssm_point_data_name is None:
        ssm_point_data_name = _default_point_data_name

    if should_read_point_data is None and isinstance(point_data_name, str):
        should_read_point_data = True
    if ssm_should_read_point_data is None and isinstance(ssm_point_data_name, str):
        ssm_should_read_point_data = True

    if isinstance(number_of_components, str):
        if number_of_components == "_none_":
            number_of_components = None
        else:
            try:
                number_of_components = int(number_of_components)
            except ValueError:
                number_of_components = None

    if structure_id_list is None:
        structure_id_list = [0, 1, 2]  # for example

    if structure_name_list is None:
        structure_name_list = _strcture_name_list_without_skin
    if structure_name_list is None:
        structure_name_list = _name2id

    _, input_basename, _ = composite_path(filename=input_file)
    mesh_reader = MeshReader(
        in_file=input_file,
        point_data_name=point_data_name,
        point_data_name_list=point_data_name_list,
        should_read_point_data=should_read_point_data,
        should_read_point_data_list=should_read_point_data_list)
    mesh = mesh_reader.run()

    output_dirname, output_basename, _ = composite_path(filename=output_file)
    output_input_file = os.path.join(output_dirname, output_basename + "__" + input_basename + ".vtk")
    output_input_disassembled_file = os.path.join(output_dirname, output_basename + "__" + input_basename + "__disassembled.vtk")

    MeshVTKWriter(mesh=mesh, out_file=output_input_file).run()

    mesh_disassembler = MeshDisassembler(
        mesh=mesh,
        structure_id_list=structure_id_list,
        structure_name_list=structure_name_list,
        name2id=name2id)
    disassembled_mesh = mesh_disassembler.run()
    MeshVTKWriter(mesh=disassembled_mesh, out_file=output_input_disassembled_file).run()

    ### SSM ###
    ssm_mean_mesh = None
    if check_file(ssm_mean_mesh_file):
        ssm_mean_mesh_reader = MeshReader(
            in_file=ssm_mean_mesh_file,
            point_data_name=ssm_point_data_name,
            point_data_name_list=ssm_point_data_name_list,
            should_read_point_data=ssm_should_read_point_data,
            should_read_point_data_list=ssm_should_read_point_data_list)
        ssm_mean_mesh = ssm_mean_mesh_reader.run()

    components = None
    if check_file(ssm_components_file):
        components = read_matrix_sitk(in_file=ssm_components_file, dtype=None)

    variances = None
    if check_file(ssm_variances_file):
        variances = read_matrix_sitk(in_file=ssm_variances_file, dtype=None)

    shape_parameter = None
    if check_file(shape_parameter_file):
        shape_parameter = read_matrix_sitk(in_file=shape_parameter_file, dtype=None)

    if shape_parameter is not None and isinstance(shape_parameter, (list, tuple)):
        shape_parameter = np.array(shape_parameter)

    if shape_parameter is not None and isinstance(shape_parameter, (np.ndarray, list, tuple)):
        number_of_components = len(shape_parameter)

    if number_of_components is not None:
        if components is not None:
            components = components[0:number_of_components]
        if variances is not None:
            variances = variances[0:number_of_components]
        if shape_parameter is not None:
            shape_parameter = shape_parameter[0:number_of_components]

    ssm_calculator = None
    if ssm_mean_mesh is not None:
        ssm_calculator = SSM(
            mean_points=ssm_mean_mesh.points,
            components=components,
            standard_deviations=np.sqrt(variances),
            mean_points_flatten=None,
            number_of_batch=None,
            logger=logger)

        output_ssm_mean_mesh_file = os.path.join(output_dirname, output_basename + "__ssm.vtk")
        output_ssm_mean_mesh_disassembled_file = os.path.join(output_dirname, output_basename + "__ssm__disassembled.vtk")

        MeshVTKWriter(mesh=ssm_mean_mesh, out_file=output_ssm_mean_mesh_file).run()

        ssm_mean_mesh_disassembler = MeshDisassembler(
            mesh=ssm_mean_mesh,
            structure_id_list=structure_id_list,
            structure_name_list=structure_name_list,
            name2id=name2id)
        disassembled_ssm_mean_mesh = ssm_mean_mesh_disassembler.run()
        MeshVTKWriter(mesh=disassembled_ssm_mean_mesh, out_file=output_ssm_mean_mesh_disassembled_file).run()

    if ssm_calculator is not None:
        output_ssm2input_deformed_txt_file = os.path.join(output_dirname, output_basename + "__" + input_basename + "__deformed.txt")
        output_ssm2input_deformed_file = os.path.join(output_dirname, output_basename + "__" + input_basename + "__deformed.vtk")
        output_ssm2input_deformed_disassembled_file = os.path.join(output_dirname, output_basename + "__" + input_basename + "__deformed__disassembled.vtk")

        shape_parameter_of_mesh = ssm_calculator.calculate_shape_parameter(points=mesh.points)
        deformed_ssm2input_mesh, disassembled_deformed_ssm2input_mesh = generate_deformed_meshes(
            shape_parameter=shape_parameter_of_mesh,
            ssm_mean_mesh=ssm_mean_mesh,
            ssm_calculator=ssm_calculator,
            structure_id_list=structure_id_list,
            structure_name_list=structure_name_list,
            name2id=name2id)

        write_matrix(out_file=output_ssm2input_deformed_txt_file, x=shape_parameter_of_mesh)
        if deformed_ssm2input_mesh is not None:
            MeshVTKWriter(mesh=deformed_ssm2input_mesh, out_file=output_ssm2input_deformed_file).run()
        if disassembled_deformed_ssm2input_mesh is not None:
            MeshVTKWriter(mesh=disassembled_deformed_ssm2input_mesh, out_file=output_ssm2input_deformed_disassembled_file).run()

    if shape_parameter is not None:
        output_ssm_deformed_file = os.path.join(output_dirname, output_basename + "__deformed.vtk")
        output_ssm_deformed_disassembled_file = os.path.join(output_dirname, output_basename + "__deformed__disassembled.vtk")

        deformed_ssm_mesh, disassembled_deformed_ssm_mesh = generate_deformed_meshes(
            shape_parameter=shape_parameter,
            ssm_mean_mesh=ssm_mean_mesh,
            ssm_calculator=ssm_calculator,
            structure_id_list=structure_id_list,
            structure_name_list=structure_name_list,
            name2id=name2id)

        if deformed_ssm_mesh is not None:
            MeshVTKWriter(mesh=deformed_ssm_mesh, out_file=output_ssm_deformed_file).run()
        if disassembled_deformed_ssm_mesh is not None:
            MeshVTKWriter(mesh=disassembled_deformed_ssm_mesh, out_file=output_ssm_deformed_disassembled_file).run()


if __name__ == '__main__':
    main(sys_argv=sys.argv)
