#
import sys
import os
import numpy as np
from collections import OrderedDict
#
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
#
dsa = None
from vtk.numpy_interface import dataset_adapter as dsa
#
torch = None
try:
    import torch
except:
    if True:
        print(("exception : import torch"))
#
sitk = None
try:
    import SimpleITK as sitk
except:
    if True:
        print(("exception : import SimpleITK as sitk"))
#
from resource import _name2id
from resource import _strcture_name_list_without_skin
from resource import _strcture_name_list
#
_default_point_data_name=("label")
#
def composite_path(filename):
    dirname, basefile = os.path.split(filename)
    basename, ext = os.path.splitext(basefile)
    return dirname, basename, ext
#
def check_file(in_file):
    if in_file is not None:
        if isinstance(in_file,(str,)):
            if os.path.exists(in_file):
                if os.path.isfile(in_file):
                    return True
    #
    return False
#
def check_folder(in_folder):
    if in_folder is not None:
        if isinstance(in_folder,(str,)):
            if os.path.exists(in_folder):
                if os.path.isdir(in_folder):
                    return True
    #
    return False
#
def read_matrix(in_file,dtype=None):
    if dtype is None:
        dtype = ("float")
    matrix=None
    if check_file(in_file):
        in_folder,in_fname,in_ext=composite_path(filename=in_file)
        if in_ext in ((".txt"),(".dat")):
            matrix = np.loadtxt(fname=in_file, dtype=dtype, comments=("#"), delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
        elif (".npy") == in_ext:
            matrix = np.load(in_file)
    #
    return matrix
    #
#
def read_matrix_sitk(in_file,dtype=None):
    if dtype is None:
        dtype = ("float")
    matrix = None
    if check_file(in_file):
        in_folder,in_fname,in_ext=composite_path(filename=in_file)
        if in_ext in ((".txt"),(".dat")):
            matrix = np.loadtxt(fname=in_file, dtype=dtype, comments=("#"), delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
        elif (".npy") == in_ext:
            matrix = np.load(in_file)
        elif in_ext in ((".mha"),(".mhd"),(".gz"),(".nii.gz"),):
            image_sitk = sitk.ReadImage(in_file)
            matrix = sitk.GetArrayFromImage(image_sitk)
    #
    return matrix
    #
#
def write_matrix(out_file,x,fmt=None,do_calc_fmt_in=None):
    #
    out_folder = None
    out_fname = None
    out_ext = None
    if out_file is not None:
        if isinstance(out_file,(str,)):
            out_folder,out_fname,out_ext = composite_path(out_file)
    #
    if isinstance(x,(dict,)):
        out_npz_file = os.path.join( out_folder,out_fname + (".npz") )
        #np.savez(out_npz_file,x)#do not
        np.savez(out_npz_file,**x)
    else:
        if do_calc_fmt_in is not None:
            if do_calc_fmt_in:
                x_dtype = x.dtype
                if x_dtype in (np.float16,np.float32,np.float64,np.float8):
                    fmt=("%lf")
                else:
                    fmt=("%d")
        #
        if fmt is None:
            fmt=("%lf")
        #
        if (".txt") == out_ext:
            np.savetxt(fname=out_file,X=x,fmt=fmt,delimiter=("\t"),newline=("\n"),header=(""),footer=(""),comments=("#"))
        elif (".npy") == out_ext:
            np.save(out_file,x)
        #
    #
#
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
    #
    #
#
class MeshReader():
    def __init__(self, in_file,point_data_name=None, point_data_name_list=None,should_read_point_data=None,should_read_point_data_list=None):
        self._in_file = in_file
        self._point_data_name = point_data_name
        self._point_data_name_list = point_data_name_list
        self._should_read_point_data = should_read_point_data
        self._should_read_point_data_list = should_read_point_data_list
        self._mesh = None
    #
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
    
    def run(self,):
        self._mesh = self.read( in_file=self.in_file, point_data_name=self.point_data_name, point_data_name_list=self.point_data_name_list )
        return self.mesh

    def read(self,in_file,point_data_name, point_data_name_list=None):
        vtk_polydata = self.read_polygon(in_file=in_file)
        if vtk_polydata is None:
            return None
        points = self.get_points(vtk_polydata=vtk_polydata)
        faces = self.get_faces(vtk_polydata=vtk_polydata)
        vtk_polydata_dsa = None
        if ( (point_data_name is not None) and isinstance(point_data_name,(str,)) ) or ( (point_data_name_list is not None) and isinstance(point_data_name_list,(list,tuple,dict,str)) ):
            vtk_polydata_dsa = self.create_vtk_polydata_dsa(vtk_polydata=vtk_polydata)
        #
        point_data = None
        point_data_dict = None
        if (vtk_polydata_dsa is not None):
            if (self.should_read_point_data is not None) and self.should_read_point_data:
                if ( (point_data_name is not None) and isinstance(point_data_name,(str,)) ):
                    point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=point_data_name)
            #
            if (self.should_read_point_data_list is not None) and self.should_read_point_data_list:
                if point_data_name_list is not None:
                    if isinstance(point_data_name_list,(str,)):
                        point_data_raw = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=point_data_name_list)
                        if point_data_raw is not None:
                            if point_data_dict is None:
                                point_data_dict = OrderedDict({})
                            point_data_dict[point_data_name_list] = point_data_raw
                    elif isinstance(point_data_name_list,(list,tuple,)):
                        point_data_dict = self.read_point_data_list(vtk_polydata_dsa=vtk_polydata_dsa, name_list=point_data_name_list)
                    elif isinstance(point_data_name_list,(dict,)):
                        point_data_dict = self.read_point_data_dict(vtk_polydata_dsa=vtk_polydata_dsa, name_dict=point_data_name_list)
        
        mesh = Mesh(
            points=points, 
            faces=faces, 
            point_data=point_data, 
            point_data_dict=point_data_dict,
            point_data_name=point_data_name
        )

        return mesh
    #
    def read_polygon(self,in_file):
        in_dirname, in_basename, in_ext = self.composite_path(filename=in_file)
        #
        reader = None
        #
        scalar = None
        if in_ext ==('.stl'):
            reader = vtk.vtkSTLReader()
        elif in_ext == ('.ply'):
            reader = vtk.vtkPLYReader()
        elif in_ext == ('.obj'):
            reader = vtk.vtkOBJReader()
        elif in_ext == ('.vtp'):
            reader = vtk.vtkXMLPolyDataReader()
            scalar = True
        elif in_ext == ('.vtk'):
            reader = vtk.vtkPolyDataReader()
        else:
            return None
        #
        if reader is not None:
            reader.SetFileName(in_file)
        #
        if reader is not None:
            reader.Update()
        #
        polydata = None
        if reader is not None:
            polydata = reader.GetOutput()
        #
        return polydata
    #
    def get_points(self,vtk_polydata):
        vtk_points = vtk_polydata.GetPoints()
        vtkarray_points = vtk_points.GetData()
        points = vtk_to_numpy(vtkarray_points)
        return points
    #
    def get_faces(self,vtk_polydata):
        vtk_cells = vtk_polydata.GetPolys()
        n_cells = vtk_cells.GetNumberOfCells()
        vtkarray_cells = vtk_cells.GetData()
        n_cols = vtkarray_cells.GetNumberOfValues()//n_cells
        cells = vtk_to_numpy(vtkarray_cells)
        cells = cells.reshape((-1,n_cols))
        cells_original = cells
        cells = cells[:,1:]
        return cells
    #
    def create_vtk_polydata_dsa(self,vtk_polydata):
        return dsa.WrapDataObject(vtk_polydata)
    #
    def get_point_data(self,vtk_polydata_dsa, name):
        return vtk_polydata_dsa.PointData[name]
    #
    def get_field_data(self,vtk_polydata_dsa, name):
        return vtk_polydata_dsa.FieldData[name]
    #
    def read_point_data(self,vtk_polydata_dsa, name):
        point_data = None
        if self.check_file(name):
            point_data = self.read_matrix(in_file=name,dtype=None)
        else:
            point_data = self.get_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=name)
        return point_data
    #
    def read_point_data_list(self,vtk_polydata_dsa, name_list):
        point_data_dict = None
        for i in range(0,name_list):
            name = name_list[i]
            if name is None:
                continue
            point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=name)
            if point_data is None:
                continue
            if point_data_dict is None:
                point_data_dict = OrderedDict({})
            point_data_dict[name] = point_data
        return point_data_dict
    #
    def read_point_data_dict(self,vtk_polydata_dsa, name_dict):
        point_data_dict = None
        for k in name_dict.keys():
            if k is None:
                continue
            v = name_dict.get(k)
            if True:
                if v is None:
                    v=k
            if v is None:
                continue
            point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=v)
            if point_data is None:
                point_data = self.read_point_data(vtk_polydata_dsa=vtk_polydata_dsa, name=k)
            if point_data is None:
                continue
            if point_data_dict is None:
                point_data_dict = OrderedDict({})
            point_data_dict[k] = point_data
        return point_data_dict
    #
    def composite_path(self,filename):
        return composite_path(filename=filename)
    #
    def check_file(self,in_file):
        return check_file(in_file=in_file)
    #
    def check_folder(self,in_folder):
        return check_folder(in_folder=in_folder)
    #
    def read_matrix(self,in_file,dtype=None):
        #return read_matrix(in_file=in_file,dtype=dtype)
        return read_matrix_sitk(in_file=in_file,dtype=dtype)
    #
#
class MeshDisassembler():
    def __init__(self, mesh, structure_id_list, structure_name_list=None, name2id=None):
        self._mesh = mesh
        self._structure_id_list = structure_id_list
        self._structure_name_list = structure_name_list
        self._name2id = name2id
        self._disassembled_mesh = None

        if self._structure_id_list is None:
            if ( (structure_name_list is not None) and isinstance(structure_name_list,(list,tuple)) ) and ( (name2id is not None) and isinstance(name2id,(dict,)) ):
                self._structure_id_list = self.make_structure_id_list(structure_name_list=structure_name_list, name2id=name2id)

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
        self._disassembled_mesh = self.process( mesh=self.mesh, structure_id_list=self.structure_id_list )
        return self.disassembled_mesh
    
    def process(self,mesh,structure_id_list):
        if structure_id_list is None:
            return None
        #
        if not isinstance(structure_id_list,(np.ndarray,tuple,list,dict,)):
            return None
        
        disassembled_mesh = None
        if isinstance(structure_id_list,(dict,)):
            name2id_table = structure_id_list
            value_set = set(name2id_table.values())
            value_list = list(value_set)
            structure_id_list = value_list
            if True:
                structure_id_list = sorted(value_list)
            disassembled_mesh = self.assemble(mesh=mesh,structure_id_list=structure_id_list)
        elif isinstance(structure_id_list,(np.ndarray,tuple,list,)):
            disassembled_mesh = self.assemble(mesh=mesh,structure_id_list=structure_id_list)
        
        return disassembled_mesh
    
    def disassemble(self,structure_id):
        self._disassembled_mesh = self.decompose(mesh=self.mesh,structure_id=structure_id)
    
    def assemble(self,mesh,structure_id_list):
        if mesh is None:
            return None
        if structure_id_list is None:
            return None
        len_of_structure_id_list = len(structure_id_list)
        mesh_points_max = 0
        points_list = None
        faces_list = None
        point_data_list = None
        point_data_dict_list = None
        for i in range(0,len_of_structure_id_list):
            structure_id = structure_id_list[i]
            decomposed_mesh = self.decompose(mesh=mesh,structure_id=structure_id)
            #
            len_of_decomposed_mesh_points = len(decomposed_mesh.points)
            #
            decomposed_mesh_faces_updated = decomposed_mesh.faces + mesh_points_max
            #
            if points_list is None:
                points_list = []
            points_list.append( decomposed_mesh.points )
            if faces_list is None:
                faces_list = []
            faces_list.append( decomposed_mesh_faces_updated )
            if point_data_list is None:
                point_data_list = []
            point_data_list.append( decomposed_mesh.point_data )
            #
            if (decomposed_mesh.point_data_dict is not None) and isinstance(decomposed_mesh.point_data_dict,(dict,)):
                for k in mesh.point_data_dict.keys():
                    if k is None:
                        continue
                    v = decomposed_mesh.point_data_dict.get(k)
                    if v is None:
                        continue
                    if point_data_dict_list is None:
                        point_data_dict_list = OrderedDict({})
                    u = point_data_dict_list.get(k)
                    if u is None:
                        point_data_dict_list[k] = []
                    point_data_dict_list[k].append(v)
            mesh_points_max = mesh_points_max + len_of_decomposed_mesh_points
        
        points_concat = None
        if points_list is not None:
            if len(points_list) > 0:
                points_concat = np.concatenate(points_list,axis=0)
        
        faces_concat = None
        if faces_list is not None:
            if len(faces_list) > 0:
                faces_concat = np.concatenate(faces_list,axis=0)
        
        point_data_concat = None
        if point_data_list is not None:
            if len(point_data_list) > 0:
                point_data_concat = np.concatenate(point_data_list,axis=0)
        #
        point_data_dict_concat = None
        #
        if (point_data_dict_list is not None) and isinstance(point_data_dict_list,(dict,)):
            for k in point_data_dict_list.keys():
                if k is None:
                    continue
                v = point_data_dict_list.get(k)
                if v is None:
                    continue
                if len(v) > 0:
                    if point_data_dict_concat is None:
                        point_data_dict_concat = OrderedDict({})
                    point_data_dict_concat[k] = np.concatenate(v,axis=0)
        #
        concat_mesh = Mesh(
            points=points_concat, 
            faces=faces_concat, 
            point_data=point_data_concat, 
            point_data_dict=point_data_dict_concat,
            point_data_name=mesh.point_data_name
        )
        return concat_mesh
        #
            
    def decompose(self,mesh,structure_id):
        if mesh is None:
            return None
        if mesh.point_data is None:
            return None
        #
        area = (mesh.point_data == structure_id)
        #
        count_of_area = np.count_nonzero(area)
        #
        if True:
            if count_of_area <= 0:
                return None
        #
        indices = np.where(area)
        #
        pts = mesh.points[area]
        pd = None
        if mesh.point_data is not None:
            pd = mesh.point_data[area]
        #
        indices_dict = self.convert_list_to_dict_tuple0(indices, do_reverse_key_value=True)#
        #
        fcs = self.extract_faces(valid_vertex_id_table=indices_dict, faces=mesh.faces)
        #
        pd_dict=None
        if (mesh.point_data_dict is not None) and isinstance(mesh.point_data_dict,(dict,)):
            for k in mesh.point_data_dict.keys():
                if k is None:
                    continue
                v = mesh.point_data_dict.get(k)
                if v is None:
                    continue
                vpd = v[area]
                if vpd is None:
                    continue
                #
                if pd_dict is None:
                    pd_dict=OrderedDict({})
                pd_dict[k] = vpd
        #
        decomposed_mesh = Mesh(
            points=pts, 
            faces=fcs, 
            point_data=pd, 
            point_data_dict=pd_dict,
            point_data_name=mesh.point_data_name
        )
        return decomposed_mesh
    
    def convert_list_to_dict(self, list_like_object, do_reverse_key_value=None):
        if do_reverse_key_value is None:
            #do_reverse_key_value=False
            do_reverse_key_value=True
        dict_object = OrderedDict({})
        len_of_list_like_object = len(list_like_object)
        for i in range(0,len_of_list_like_object):
            value = list_like_object[i]
            if do_reverse_key_value:
                dict_object[value] = i
            else:
                dict_object[i] = value
        #
        return dict_object
    #
    def convert_list_to_dict_tuple0(self, tuple_list_like_object, do_reverse_key_value=None):
        return self.convert_list_to_dict(list_like_object=tuple_list_like_object[0], do_reverse_key_value=do_reverse_key_value)
    #
    def extract_faces(self, valid_vertex_id_table, faces):#
        len_of_faces = len(faces)
        new_face_list = []
        for i in range(0,len_of_faces):
            face = faces[i]
            face_0 = valid_vertex_id_table.get(face[0])
            face_1 = valid_vertex_id_table.get(face[1])
            face_2 = valid_vertex_id_table.get(face[2])
            #
            if (face_0 is not None) and (face_1 is not None) and (face_2 is not None):
                new_fc = (face_0,face_1,face_2)
                new_face = np.array(new_fc)
                new_face_list.append(new_face)
            #
        #
        return np.array(new_face_list)
        #
    #
    def make_structure_id_list(self,structure_name_list, name2id):
        structure_id_list=None
        for i in range(0,len(structure_name_list)):
            structure_name = structure_name_list[i]
            if structure_name is None:
                continue
            structure_id = name2id.get(structure_name)
            if structure_id is None:
                continue
            if structure_id_list is None:
                structure_id_list.append( structure_id )
        #
        return structure_id_list
        #
#
class MeshVTKWriter():
    def __init__(self, mesh, out_file):
        self._mesh = mesh
        self._out_file = out_file
    #
    @property
    def mesh(self):
        return self._mesh

    @property
    def out_file(self):
        return self._out_file
    
    def write_vtk_polygon_pts(self, ff,points):
        if points is None:
            return
        len_of_points = len(points)
        #
        ff.write(("POINTS {} float\n").format(len_of_points))
        #
        for i in range(0,len_of_points):
            point = points[i]
            ff.write(("{} {} {} \n").format(point[0],point[1],point[2]))
        #
    #
    def write_vtk_polygon_faces(self,ff,faces,mesh_ndim=3):
        if faces is None:
            return
        faces_shape = faces.shape
        len_of_faces = faces_shape[0]
        dim_of_content = (faces_shape[0] * faces_shape[1]) + len_of_faces
        #
        ff.write(("POLYGONS {} {}\n").format(len_of_faces,dim_of_content))
        #
        for i in range(0,len_of_faces):
            face = faces[i]
            ff.write(("{} {} {} {} \n").format(mesh_ndim,face[0],face[1],face[2]))
        #
    #
    def write_vtk_polygon_scalar_vector(self,ff,name,scalar_vector):
        if scalar_vector is None:
            return
        #
        if name is None:
            name = ("scalar")
        #
        len_of_scalar_vector = len(scalar_vector)
        scalar_vector_dtype = scalar_vector.dtype
        #
        #
        if scalar_vector_dtype == np.float64:
            ff.write(("{} 1 {} double\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.float32:
            ff.write(("{} 1 {} float\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.int32:
            ff.write(("{} 1 {} int\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.uint32:
            ff.write(("{} 1 {} unsigned_int\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.int8:
            ff.write(("{} 1 {} char\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.uint8:
            ff.write(("{} 1 {}  unsigned_char\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.int16:
            ff.write(("{} 1 {} short\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.uint16:
            ff.write(("{} 1 {}  unsigned_short\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.int64:
            if False:
                if True:
                    scalar_vector = scalar_vector.astype(np.int32)
                #
                ff.write(("{} 1 {} int\n").format(name,len_of_scalar_vector))
            else:
                ff.write(("{} 1 {} long\n").format(name,len_of_scalar_vector))
        elif scalar_vector_dtype == np.uint64:
            ff.write(("{} 1 {} unsigned_long\n").format(name,len_of_scalar_vector))
        else:
            ff.write(("{} 1 {} double\n").format(name,len_of_scalar_vector))
        #
        for i in range(0,len_of_scalar_vector):
            scalar = scalar_vector[i]
            try:
                len_of_scalar = len(scalar)
                text = ("")
                for j in range(0,len_of_scalar):
                    text += str(scalar[j]) + (" ")
                text += ("\n")
                ff.write(text)
            except:
                ff.write(("{} \n").format(scalar))
            #
        #
    #
    def _write_vtk_polygon_file(self,filename,points,faces,scalar_vector=None,scalar_vector_name=None,scalar_vector_dict=None,do_cast_points=None,mesh_ndim=None,write_mode=None):
        if filename is None:
            return False
        #
        if points is None:
            return False
        #
        if faces is None:
            return False
        #
        if do_cast_points is None:
            do_cast_points=True
        if mesh_ndim is None:
            if points is None:
                mesh_ndim=3
            else:
                mesh_ndim = points.shape[1]
        if write_mode is None:
            write_mode=("w")
        #
        if True:
        #if False:
            if do_cast_points:
                if points is not None:
                    points = points.astype(np.float32)
                #
            #
        #
        if mesh_ndim is None:
            mesh_ndim = points.shape[1]#3
        #
        len_of_points = len(points)
        #
        with open(filename, write_mode) as ff:
            ff.write(("# vtk DataFile Version 4.2\n"))
            ff.write(("vtk output\n"))
            #
            if ("wb") == write_mode:
                ff.write(("BINARY\n"))
            else:
                ff.write(("ASCII\n"))
            #
            ff.write(("DATASET POLYDATA\n"))
            self.write_vtk_polygon_pts(ff=ff,points=points)
            self.write_vtk_polygon_faces(ff=ff,faces=faces,mesh_ndim=mesh_ndim)
            #
            if scalar_vector is not None or ( scalar_vector_dict is not None and isinstance(scalar_vector_dict,(dict,)) ):
                ff.write(("POINT_DATA {}\n").format(len_of_points))
                counter = 0
                len_of_scalar_vector = None
                if scalar_vector is not None:
                    len_of_scalar_vector = len(scalar_vector)
                    counter = counter + 1
                #
                scalar_vector_dict_key_list = None
                len_of_scalar_vector_dict_key_list = None
                if scalar_vector_dict is not None:
                    if isinstance(scalar_vector_dict,(dict,)):
                        scalar_vector_dict_key_list = list(scalar_vector_dict.keys())
                        len_of_scalar_vector_dict_key_list = len(scalar_vector_dict_key_list)
                        counter = counter + len_of_scalar_vector_dict_key_list
                #
                ff.write(("FIELD FieldData {}\n").format(counter))
                #
                if len_of_scalar_vector is not None:
                    self.write_vtk_polygon_scalar_vector(ff=ff,name=scalar_vector_name,scalar_vector=scalar_vector)
                #
                if len_of_scalar_vector_dict_key_list is not None:
                    for i in range(0,len_of_scalar_vector_dict_key_list):
                        scalar_vector_dict_key = scalar_vector_dict_key_list[i]
                        scalar_vector_dict_value = scalar_vector_dict.get(scalar_vector_dict_key)
                        if scalar_vector_dict_value is not None:
                            len_of_scalar_vector_dict_value = len(scalar_vector_dict_value)
                            if len_of_points == len_of_scalar_vector_dict_value:
                                self.write_vtk_polygon_scalar_vector(ff=ff,name=scalar_vector_dict_key,scalar_vector=scalar_vector_dict_value)
                #
                
            #
            if False:
                if scalar_vector is not None:
                    len_of_scalar_vector = len(scalar_vector)
                    #len_of_points = len(points)
                    if len_of_points == len_of_scalar_vector:
                        self.write_vtk_polygon_scalar_vector(ff=ff,name=scalar_vector_name,scalar_vector=scalar_vector)
                #
                if scalar_vector_dict is not None:
                    if isinstance(scalar_vector_dict,(dict,)):
                        scalar_vector_dict_key_list = list(scalar_vector_dict.keys())
                        len_of_scalar_vector_dict_key_list = len(scalar_vector_dict_key_list)
                        for i in range(0,len_of_scalar_vector_dict_key_list):
                            scalar_vector_dict_key = scalar_vector_dict_key_list[i]
                            scalar_vector_dict_value = scalar_vector_dict.get(scalar_vector_dict_key)
                            if scalar_vector_dict_value is not None:
                                len_of_scalar_vector_dict_value = len(scalar_vector_dict_value)
                                if len_of_points == len_of_scalar_vector_dict_value:
                                    self.write_vtk_polygon_scalar_vector(ff=ff,name=scalar_vector_dict_key,scalar_vector=scalar_vector_dict_value)
                #
            #
        #
        return True
        #
    #
    def write_vtk_polygon_file(self,do_cast_points=None,mesh_ndim=None,write_mode=None):
        out_file = self.out_file
        out_dirname, out_basename, out_ext = self.composite_path(filename=out_file)
        os.makedirs(out_dirname,exist_ok=True)
        return self._write_vtk_polygon_file(filename=out_file,points=self.mesh.points,faces=self.mesh.faces,scalar_vector=self.mesh.point_data,scalar_vector_name=self.mesh.point_data_name,scalar_vector_dict=self.mesh.point_data_dict,do_cast_points=do_cast_points,mesh_ndim=mesh_ndim,write_mode=write_mode)
    #
    def run(self,):
        do_cast_points=None
        mesh_ndim=None
        write_mode=None
        return self.write_vtk_polygon_file(do_cast_points=do_cast_points,mesh_ndim=mesh_ndim,write_mode=write_mode)
    #
    def composite_path(self,filename):
        return composite_path(filename=filename)
    #
#
#
class SSM():
    def __init__(self,
        mean_points,
        components,
        standard_deviations,
        mean_points_flatten=None,
        number_of_batch=None,
        logger=None,
    ):
        self._mean_points = mean_points
        if mean_points_flatten is None:
            if mean_points is not None:
                #mean_points_flatten = torch.flatten(mean_points)
                #mean_points_flatten = mean_points.ravel()
                #mean_points_flatten = mean_points.flatten()
                if isinstance(mean_points,(np.ndarray,)):
                    mean_points_flatten = mean_points.ravel()
                else:
                    mean_points_flatten = torch.flatten(mean_points)
        self._mean_points_flatten = mean_points_flatten
        self._components = components
        self._standard_deviations = standard_deviations
        #
        if number_of_batch is None:
            number_of_batch=1
        self._number_of_batch = number_of_batch
        #
        mean_points_batch = None
        #
        self._mean_points_batch = None
        self._mean_points_batch = mean_points_batch
        self._mean_points_batch_flatten = None
        #
        self._logger = logger
    #
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
    #
    @property
    def mean_points_flatten(self):
        return self._mean_points_flatten
    #
    def _deform_term(self,shape_parameter,standard_deviations=None):
        if standard_deviations is None:
            standard_deviations = self.standard_deviations
        shapre_score = standard_deviations*shape_parameter
        deform_flatten = None
        if isinstance(self.components,(np.ndarray,)):
            deform_flatten = np.matmul(shapre_score,self.components)
        else:
            deform_flatten = torch.matmul(shapre_score,self.components)
        return deform_flatten
    #
    def _deform(self,shape_parameter,standard_deviations=None,mean_points_flatten=None):
        if mean_points_flatten is None:
            mean_points_flatten = self.mean_points_flatten
        deform_term = self._deform_term(shape_parameter=shape_parameter,standard_deviations=standard_deviations)
        deformed_verts_flatten = mean_points_flatten + deform_term
        return deformed_verts_flatten
    #
    def deform(self,shape_parameter):
        if shape_parameter is None:
            return self.mean_points
        standard_deviations = self.standard_deviations
        mean_points_flatten = self.mean_points_flatten
        deformed_verts_flatten=self._deform(shape_parameter=shape_parameter,standard_deviations=standard_deviations,mean_points_flatten=mean_points_flatten)
        return deformed_verts_flatten.reshape(-1,self.mean_points.shape[1])
    #
    def deform_batch(self,shape_parameter):
        if shape_parameter is None:
            return self._mean_points_batch
        standard_deviations = None
        if isinstance(self.standard_deviations,(np.ndarray,)):
            standard_deviations = np.expand_dims(self.standard_deviations, axis=0)
        else:
            standard_deviations = torch.unsqueeze(self.standard_deviations, dim=0)
        #
        mean_points_flatten = self.mean_points_flatten
        mean_points_flatten_true = mean_points_flatten
        deformed_verts_flatten=self._deform(shape_parameter=shape_parameter,standard_deviations=standard_deviations,mean_points_flatten=mean_points_flatten_true)
        return deformed_verts_flatten.reshape(self._number_of_batch,-1,self.mean_points.shape[1])
    #
    def calculate_shape_parameter(self,points):
        points_flatten = points.ravel()
        diff = points_flatten - self.mean_points_flatten
        score = np.matmul(self.components, diff)
        normalized_score = score / self.standard_deviations
        return normalized_score
    #
#
#
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
    #
    if sys_argv is not None:
        if isinstance(sys_argv,(list,tuple)):
            len_of_sys_argv = len(sys_argv)
            if len_of_sys_argv>0:
                if mode is None:
                    mode = 0
                if 0 == mode:
                    if len_of_sys_argv==2:
                        input_file = str(sys_argv[1])
                    elif len_of_sys_argv==3:
                        input_file = str(sys_argv[1])
                        output_file = str(sys_argv[2])
                    elif len_of_sys_argv==4:
                        input_file = str(sys_argv[1])
                        output_file = str(sys_argv[2])
                        point_data_name = str(sys_argv[3])
                    elif len_of_sys_argv==8:
                        input_file = str(sys_argv[1])
                        output_file = str(sys_argv[2])
                        point_data_name = str(sys_argv[3])
                        ssm_mean_mesh_file = str(sys_argv[4])
                        ssm_components_file = str(sys_argv[5])
                        ssm_variances_file = str(sys_argv[6])
                        number_of_components = str(sys_argv[7])
                    elif len_of_sys_argv==9:
                        input_file = str(sys_argv[1])
                        output_file = str(sys_argv[2])
                        point_data_name = str(sys_argv[3])
                        ssm_mean_mesh_file = str(sys_argv[4])
                        ssm_components_file = str(sys_argv[5])
                        ssm_variances_file = str(sys_argv[6])
                        number_of_components = str(sys_argv[7])
                        shape_parameter_file = str(sys_argv[8])
                    else:
                        print ("1:input_file(str)[.vtk]")
                        print ("2:output_file(str)[.vtk]")
                        print ("3:point_data_name(str)")
                        print ("4:ssm_mean_mesh_file(str)[.vtk]")
                        print ("5:ssm_components_file(str)[.txt,.mhd/.mha]")
                        print ("6:ssm_variances_file(str)[.txt]")
                        print ("7:number_of_components(str,int)","_none_")
                        print ("8:shape_parameter_file(str)[.txt]")
                        return
    #
    if True:
        print(("input_file"),input_file)
        print(("output_file"),output_file)
        print(("point_data_name"),point_data_name)
        print(("ssm_mean_mesh_file"),ssm_mean_mesh_file)
        print(("ssm_components_file"),ssm_components_file)
        print(("ssm_variances_file"),ssm_variances_file)
        print(("number_of_components"),number_of_components)
        print(("shape_parameter_file"),shape_parameter_file)
    #
    if point_data_name is None:
        point_data_name=_default_point_data_name
    #
    if ssm_point_data_name is None:
        ssm_point_data_name=_default_point_data_name
    #
    if should_read_point_data is None:
        if ( (point_data_name is not None) and isinstance(point_data_name,(str,)) ):
            should_read_point_data=True
    #
    if ssm_should_read_point_data is None:
        if ( (ssm_point_data_name is not None) and isinstance(ssm_point_data_name,(str,)) ):
            ssm_should_read_point_data=True
    #
    if number_of_components is not None:
        if isinstance(number_of_components,(str,)):
            if ("_none_") == number_of_components:
                number_of_components=None
            else:
                try:
                    number_of_components = int(number_of_components)
                except:
                    number_of_components = None
                    pass
    #
    if True:
        if structure_id_list is None:
            structure_id_list = [0,1,2] # for example
            #structure_id_list = None
    #
    if structure_name_list is None:
        structure_name_list = _strcture_name_list_without_skin
    #
    if structure_name_list is None:
        structure_name_list = _name2id
    #
    input_dirname, input_basename, input_ext = composite_path(filename=input_file)
    mesh_reader = MeshReader(
        in_file=input_file,
        point_data_name=point_data_name, 
        point_data_name_list=point_data_name_list,
        should_read_point_data=should_read_point_data,
        should_read_point_data_list=should_read_point_data_list
    )
    mesh = mesh_reader.run()
    
    output_dirname, output_basename, output_ext = composite_path(filename=output_file)
    output_input_file = os.path.join(output_dirname,output_basename+("__")+input_basename+(".vtk"))
    output_input_disassembled_file = os.path.join(output_dirname,output_basename+("__")+input_basename+("__disassembled")+(".vtk"))

    if True:
        mesh_vtk_writer = MeshVTKWriter(mesh=mesh, out_file=output_input_file)
        mesh_vtk_writer.run()
    
    if True:
        mesh_disassembler = MeshDisassembler(
            mesh=mesh, 
            structure_id_list = structure_id_list, 
            structure_name_list=structure_name_list, 
            name2id=name2id
        )
        disassembled_mesh = mesh_disassembler.run()

        disassembled_mesh_vtk_writer = MeshVTKWriter(mesh=disassembled_mesh, out_file=output_input_disassembled_file)
        disassembled_mesh_vtk_writer.run()

    if False:
        return

    ### SSM ###
    ssm_mean_mesh = None
    if check_file(ssm_mean_mesh_file):
        ssm_mean_mesh_reader = MeshReader(
            in_file=ssm_mean_mesh_file,
            point_data_name=ssm_point_data_name, 
            point_data_name_list=ssm_point_data_name_list,
            should_read_point_data=ssm_should_read_point_data,
            should_read_point_data_list=ssm_should_read_point_data_list
        )
        ssm_mean_mesh = ssm_mean_mesh_reader.run()
    
    components = None
    if check_file(ssm_components_file):
        components = read_matrix_sitk(in_file=ssm_components_file,dtype=None)
    
    variances = None
    if check_file(ssm_variances_file):
        variances = read_matrix_sitk(in_file=ssm_variances_file,dtype=None)
    
    shape_parameter = None
    if check_file(shape_parameter_file):
        shape_parameter = read_matrix_sitk(in_file=shape_parameter_file,dtype=None)
    
    if shape_parameter is not None:
        if isinstance(shape_parameter,(list,tuple)):
            shape_parameter = np.array(shape_parameter)
    
    if shape_parameter is not None:
        if isinstance(shape_parameter,(np.ndarray,list,tuple)):
            number_of_components = len(shape_parameter) # ToDo: Error check for the number of elements is required
    #
    if number_of_components is not None:
        if components is not None:
            components = components[0:number_of_components]
        if variances is not None:
            variances = variances[0:number_of_components]
        if shape_parameter is not None:
            shape_parameter = shape_parameter[0:number_of_components]
    #
    ssm_calculator = None
    if ssm_mean_mesh is not None:
        ssm_calculator = SSM(
            mean_points=ssm_mean_mesh.points,
            components=components,
            standard_deviations=np.sqrt(variances),
            mean_points_flatten=None,
            number_of_batch=None,
            logger=logger
        )
        #
        output_ssm_mean_mesh_file = os.path.join(output_dirname,output_basename+("__")+("ssm")+(".vtk"))
        output_ssm_mean_mesh_disassembled_file = os.path.join(output_dirname,output_basename+("__")+("ssm")+("__disassembled")+(".vtk"))
        #
        if True:
            ssm_mean_mesh_vtk_writer = MeshVTKWriter(mesh=ssm_mean_mesh, out_file=output_ssm_mean_mesh_file)
            ssm_mean_mesh_vtk_writer.run()
        
        if True:
            ssm_mean_mesh_disassembler = MeshDisassembler(
                mesh=ssm_mean_mesh, 
                structure_id_list = structure_id_list, 
                structure_name_list=structure_name_list, 
                name2id=name2id
            )
            disassembled_ssm_mean_mesh = ssm_mean_mesh_disassembler.run()

            disassembled_ssm_mean_mesh_vtk_writer = MeshVTKWriter(mesh=disassembled_ssm_mean_mesh, out_file=output_ssm_mean_mesh_disassembled_file)
            disassembled_ssm_mean_mesh_vtk_writer.run()
    #
    if True:
        if ssm_calculator is not None:
            output_ssm2input_deformed_txt_file = os.path.join(output_dirname,output_basename+("__")+input_basename+("__deformed")+(".txt"))
            output_ssm2input_deformed_file = os.path.join(output_dirname,output_basename+("__")+input_basename+("__deformed")+(".vtk"))
            output_ssm2input_deformed_disassembled_file = os.path.join(output_dirname,output_basename+("__")+input_basename+("__deformed")+("__disassembled")+(".vtk"))
            #
            shape_parameter_of_mesh = ssm_calculator.calculate_shape_parameter(points=mesh.points)
            deformed_ssm2input_mesh,disassembled_deformed_ssm2input_mesh = generate_deformed_meshes(
                shape_parameter=shape_parameter_of_mesh,
                ssm_mean_mesh=ssm_mean_mesh,
                ssm_calculator=ssm_calculator,
                structure_id_list=structure_id_list, 
                structure_name_list=structure_name_list, 
                name2id=name2id
            )
            if True:
                write_matrix(out_file=output_ssm2input_deformed_txt_file,x=shape_parameter_of_mesh,fmt=None,do_calc_fmt_in=None)
            if True:
                if deformed_ssm2input_mesh is not None:
                    mesh_vtk_writer = MeshVTKWriter(mesh=deformed_ssm2input_mesh, out_file=output_ssm2input_deformed_file)
                    mesh_vtk_writer.run()
            if True:
                if disassembled_deformed_ssm2input_mesh is not None:
                    disassembled_deformed_ssm2input_mesh_vtk_writer = MeshVTKWriter(mesh=disassembled_deformed_ssm2input_mesh, out_file=output_ssm2input_deformed_disassembled_file)
                    disassembled_deformed_ssm2input_mesh_vtk_writer.run()
    #
    if True:
        if shape_parameter is not None:
            output_ssm_deformed_file = os.path.join(output_dirname,output_basename+("__deformed")+(".vtk"))
            output_ssm_deformed_disassembled_file = os.path.join(output_dirname,output_basename+("__deformed")+("__disassembled")+(".vtk"))
            #
            deformed_ssm_mesh,disassembled_deformed_ssm_mesh = generate_deformed_meshes(
                shape_parameter=shape_parameter,
                ssm_mean_mesh=ssm_mean_mesh,
                ssm_calculator=ssm_calculator,
                structure_id_list=structure_id_list, 
                structure_name_list=structure_name_list, 
                name2id=name2id
            )
            if True:
                if deformed_ssm_mesh is not None:
                    mesh_vtk_writer = MeshVTKWriter(mesh=deformed_ssm_mesh, out_file=output_ssm_deformed_file)
                    mesh_vtk_writer.run()
            if True:
                if disassembled_deformed_ssm_mesh is not None:
                    disassembled_deformed_ssm_mesh_vtk_writer = MeshVTKWriter(mesh=disassembled_deformed_ssm_mesh, out_file=output_ssm_deformed_disassembled_file)
                    disassembled_deformed_ssm_mesh_vtk_writer.run()
    #
#
def generate_deformed_meshes(
    shape_parameter,
    ssm_mean_mesh,
    ssm_calculator,
    structure_id_list, 
    structure_name_list=None, 
    name2id=None
):
    deformed_points = ssm_calculator.deform(shape_parameter=shape_parameter)
    deformed_ssm_mesh = Mesh(
        points=deformed_points, 
        faces=ssm_mean_mesh.faces, 
        point_data=ssm_mean_mesh.point_data, 
        point_data_dict=ssm_mean_mesh.point_data_dict,
        point_data_name=ssm_mean_mesh.point_data_name
    )
    #
    deformed_ssm_mesh_disassembler = MeshDisassembler(
        mesh=deformed_ssm_mesh, 
        structure_id_list = structure_id_list, 
        structure_name_list=structure_name_list, 
        name2id=name2id
    )
    disassembled_deformed_ssm_mesh = deformed_ssm_mesh_disassembler.run()
    return deformed_ssm_mesh,disassembled_deformed_ssm_mesh
#
if __name__ == '__main__':
    main(sys_argv=sys.argv)
#