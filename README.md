# LIGGGGHTS 2 GATE - Make DEM radioactive

The simple LIGGGHTS to GATE script enables the user to transform the particle data into a GATE- readable dataformat in 2 different steps.  

1. Generate a Volumetric Average of all particle positions to generate the scatter region due to other particles
2. Generate a trajectory file to use one or multiple particles as tracers

The user needs to include those files in GATE. Also a lookup-table for the different Materials in the Scatter Volume must be defined.  

## Usage
In Ipython simply run
```python3
run base.py 
l2g = Liggghts2Gate("data/",tracer_time = 100)
l2g.generate_data(  image_filename="volume.mhd",
                    traj_filename="gate.traj"
                 ) 
```

This will take all .vtk files from the `data` folder and generate the volumetric file.  
It will be saved in `image_filename` -path also the tracer trajectory will be saved in  
 `traj_filename`.  
 The variable `tracer_time` helps choosing the number of particles that will be used as  
 a tracer. Their positions will be added to a single big trajectory.
 You can also manualy choose the particle ID's. Simply choose the variable `tracer_id`  
 and pass either a single ID or a list/array.
 
If you want to have multiple tracers you need to rerun the whole script whith different  
Tracer IDs. I am sorry, Thats on our To-Do! 

## License

I am pretty sure leonard will figure something out for that!

