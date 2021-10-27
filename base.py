import numpy as np
from vtk import vtkDataSetReader
from vtk.util.numpy_support import vtk_to_numpy
from glob import glob
from natsort import natsorted
import SimpleITK as sitk
from scipy.interpolate import interp1d
class Dem2Gate():
    def __init__(self, 
                folder,
                simulation_timestep = 1e-5, 
                tracer_id = None,
                tracer_time = None,
                xlim = 50,
                ylim = 50,
                zlim = 50,
                ):
        """
        Dem2Gate Class: Use DEM outputdata and transfere it
        into a GATE readable format. Outputs two files:
        - Density image of the system
        - Tracer trajectory for Gate
        
        You still need to make a file that tells gate the density range of each pixel.  
        
        Parameters
        ----------
        folder: str
            Path to the folder containing all the vtk files which where 
            output by DEM software, e.g. LIGGGHTS.
        
        simulation_timetsep: float
            Timestep of the simulation used in DEM software.
        
        tracer_id: int,list,numpy-ndarray 
            The DEM internal particle ID of the particle used as a trajectory.
            Multiple IDs in a list/array are possible.
        
        tracer_time: float
            The time an artificial tracer should travel in the system. 
            If no tracer ID is given, the number of tracers needed to
            get to this time is calculated. The IDs used then start at 0.
            Ignored if tracer_id is not None.
            
        xlim,ylim,zlim: int
            Number of "voxels" of the output image in each dimension.
            
        """
        
        self.simulation_timestep = simulation_timestep
        self.sim_files, self. sim_times = self.get_sim_files(folder)
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        if tracer_id is None:
            if tracer_time is None:
                raise ValueError("Either a tracer_id or a tracer_time must be given")
            # get the number of tracers required to fullfill the time 
            number_of_tracers =int(tracer_time/self.sim_times[-1])+1
            self.tracer_id = np.linspace(0,number_of_tracers-1,number_of_tracers)
        else:
            if type(tracer_id) == list or type(tracer_id) == np.ndarray :
                self.tracer_id = np.asarray(tracer_id)
            else:
                self.tracer_id = np.asarray([tracer_id])
        
    def get_sim_files(self,folder):
        file_arr = [x for x in  glob(folder+"/*.vtk") if not "bound" in x]
        filenames = natsorted(file_arr,key=lambda y: y.lower())
        times = np.asarray([int(x.split(".vtk")[0].split("_")[-1]) for x in filenames]) * self.simulation_timestep
        return filenames, times

    def read_vtk(self,filename=None, timestep = None):
        if filename is None:
            if timestep is None:
                raise ValueError("No file or timestep to read given")
            filename = self.sim_files[timestep]
        reader = vtkDataSetReader()
        reader.SetFileName(filename)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        data = reader.GetOutput()
        p = data.GetPointData()
        position = vtk_to_numpy(data.GetPoints().GetData())
        radius = vtk_to_numpy(p.GetArray("radius"))
        p_id = vtk_to_numpy(p.GetArray("id"))
        return position, radius,p_id
        
    def get_dimensions(self):
        # take a random file
        file = self.sim_files[0]
        # change name to bounding box 
        file_b =  "_".join(file.split("_")[:-1])+"_boundingBox_"+file.split("_")[-1]
        with open(file_b,"r") as f:
            lines = f.readlines()
        xdim = np.asarray(lines[6].split(" ")[0:2], dtype = float)
        ydim = np.asarray(lines[8].split(" ")[0:2], dtype = float)
        zdim = np.asarray(lines[10].split(" ")[0:2], dtype = float)
        return [xdim[0],ydim[0],zdim[0]],[xdim[1],ydim[1],zdim[1]]
        
    def generate_data(self,image_filename, traj_filename, timestep = 0.0001):
        voxels = np.zeros((self.xlim, self.ylim, self.zlim))
        pos = np.asarray([[0,0,0]])
        tracer_path =	[[] for x in self.tracer_id]
        min_val, max_val = self.get_dimensions()
        for id, time in enumerate(self.sim_times):
            posi, rad, p_id = self.read_vtk(timestep = id)
            # appending tracer path
            for t_id in self.tracer_id:
                tracer_id = np.argmin(np.abs(p_id-t_id))
                tracer_path[int(t_id)].append(posi[tracer_id])
            # Make the voxel
            for position in posi:
                # find cells
                xi = int( (position[0]-min_val[0])/(max_val[0]-min_val[0]) * (self.xlim -1) )
                yi = int( (position[1]-min_val[1])/(max_val[1]-min_val[1]) * (self.ylim -1) )
                zi = int( (position[2]-min_val[2])/(max_val[2]-min_val[2]) * (self.zlim -1) )
                voxels[xi, yi, zi] += 1
        tracer_path = np.asarray(tracer_path).reshape((-1,3))
        #voxels, max_val = self.make_voxels(np.asarray(pos[1::],dtype = np.float))
        self.make_voxel_image(voxels, image_filename ,max_val, min_val)
        print(tracer_path)
        self.make_tracer_trajectory(tracer_path, traj_filename, timestep)
        return tracer_path
        
    def __repr__(self):
        string = "DEM to GATE class\n"
        return string
    
    def __str__(self):
        string = f"Currently {len(self.sim_files)} files in the Class.\n" +\
                 f"Max time is {max(self.sim_times):.1f} seconds"
        return string
    
    def __del__(self):
        pass

            
            
            
    def make_voxel_image(self, voxels,name, max_dim, min_dim):
        if not name.endswith(".mhd"):
            raise ValueError(".mhd is the only valid output type")
        img_new = sitk.GetImageFromArray(voxels)
        sitk.WriteImage(img_new, name)
        with open(name,"r") as f:
            mhd_input = f.readlines()
        mhd_input[9] = f"ElementSpacing = {(max_dim[2]-min_dim[2])/self.zlim*1000} {(max_dim[1]-min_dim[1])/self.ylim*1000} {(max_dim[0]-min_dim[0])/self.xlim*1000}\n"
        with open(name,"w") as f:
            f.writelines(mhd_input)       
                
                
    def make_tracer_trajectory(self,trajectory,filename = 'GATE-Traj.placements' , timestep = 0.0001):

        header = '''###### List of placement (translation and rotation) according to time\n\
        ###### Column 1 is Time in s (second)\n\
        ###### Column 2 is rotationAngle in degree\n\
        ###### Columns 3,4,5 are rotation axis\n\
        ###### Columns 6,7,8 are translation in mm\n\
        Time s\n\
        Rotation deg\n\
        Translation mm\n'''
        theta = 0
        axis1 = 0
        axis2 = 0
        axis3 = 0
        t = np.linspace(0,self.sim_times[-1]*len(self.tracer_id),len(trajectory[:,0]))
        x = trajectory[:,0] * 1000
        y = trajectory[:,1] * 1000
        z = trajectory[:,2] * 1000
        fx = interp1d(t, x)
        fy = interp1d(t, y)
        fz = interp1d(t, z)
        ts = np.arange(t[0], t[-1], timestep)
        xs = fx(ts)
        ys = fy(ts)
        zs = fz(ts)
        with  open(filename, 'w') as f:
            f.write(header)
            for k in range(len(ts)):
                a = str(ts[k]) + ' ' + str(theta) + ' ' + str(axis1) + ' ' + str(axis2) + ' ' + str(axis3) + \
                             ' ' + str(xs[k]) + ' ' + str(ys[k]) + ' ' + str(zs[k])+'\n'
                f.write(a)

     
