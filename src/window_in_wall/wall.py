import math
import compas

from compas.datastructures import Mesh
import itertools


class Wall(object):
    """ Wall generator

     """

    def __init__(self):
        """ constructor

        """
        self.mesh = Mesh()
        self.mesh_offset = Mesh()
        self.mesh_displaced = Mesh()
        self.gate_points = []

    def create_quad_meshes_from_dimensions(self, width=3.5, height=3.2, elsize=0.025):
        self.mesh = self.create_quad_mesh_from_dimensions(width, height, elsize)
    
    def create_quad_meshes_from_dimensions_dl(self, width=3.5, height=3.2, elsize=0.025):
        self.mesh = self.create_quad_mesh_from_dimensions(width, height, elsize)
        self.mesh_offset = self.create_quad_mesh_from_dimensions(width, height, elsize)

    def create_quad_mesh_from_dimensions(self, width=3.5, height=3.2, elsize=0.025):
        """create a regular quad mesh from height and width and edge length in the xz-plane
        """
        
        self.elsize = elsize

        self.x_size = int(width / elsize) + 1
        self.z_size = int(height / elsize) + 1

        self.width = elsize*(self.x_size-1)
        self.height = elsize*(self.z_size-1)

        mesh = Mesh()

        # create vertices
        for i in range(self.x_size):
            for j in range(self.z_size):
                x = i * self.elsize
                y = 0
                z = j * self.elsize
                glob_id = i*self.z_size+j
                mesh.add_vertex(x = x, y = 0, z = z, i=i, j=j, glob_id=glob_id, x_disp=0, z_disp=0)

        # create faces
        for i in range(self.x_size - 1):
            for j in range(self.z_size - 1):
                a = i * self.z_size + j
                b = a + self.z_size
                c = b + 1
                d = a + 1
                mesh.add_face([a, b, c, d])
        
        return mesh
    
    def sin_wave(self, amp, freq, phase, value):
        return(amp * math.sin(2.0 * math.pi * freq * value + phase))

    def undulate(self, up_amp=0.075, up_freq=2.0, up_phase=0.0,down_amp=0.075, down_freq=2.0, down_phase=0.0):
        """undulate the mesh with sin waves
        """
        for vertex in self.mesh.vertices():
            # get x coordinate from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')
            # z coordinate of the point normalized by the height of the wall
            rel_z=z_val/(self.z_size*self.elsize)
            # calculate y value of the up and down part of the wall at this x
            up_y = self.sin_wave(up_amp, up_freq, up_phase, x_val)
            down_y = self.sin_wave(down_amp, down_freq, down_phase, x_val)
            # specify a y value as a linear blending of the top and the botom of the wall
            y_val = rel_z*up_y + (1.0-rel_z)*down_y
            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=y_val)
    
    def undulate_with_displ(self, up_amp=0.075, up_freq=2.0, up_phase=0.0, down_amp=0.075, down_freq=2.0, down_phase=0.0):
        """undulate the mesh with sin waves using the displacement values
        """
        for vertex in self.mesh.vertices():
            vrt_attrs = self.mesh.vertex_attributes(vertex)
            x_new = vrt_attrs["x"] - vrt_attrs["x_disp"] #*0.99 # calc new x location
            z_new = vrt_attrs["z"] - vrt_attrs["z_disp"] #*0.99 # calc new x location

            # z coordinate of the point normalized by the height of the wall
            rel_z = z_new/(self.z_size*self.elsize)

            # calculate y value of the up and down part of the wall at this x
            up_y = self.sin_wave(up_amp, up_freq, up_phase, x_new)
            down_y = self.sin_wave(down_amp, down_freq, down_phase, x_new)

            # specify a y value as a linear blending of the top and the botom of the wall
            y_val = rel_z*up_y + (1.0-rel_z)*down_y

            # shift the y value in relation to the max amplitude:
            if len(list(self.mesh_offset.vertices())):
                layer_amp = rel_z*up_amp + (1.0-rel_z)*down_amp
                shift_val = down_amp - layer_amp
                y_val = y_val + shift_val
                self.mesh_offset.vertex_attribute(vertex, name='y', value=-y_val+2*down_amp)

            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=y_val)

    def undulate_with_displ_double(self, up_amp=0.075, up_freq=2.0, up_phase=0.0, down_amp=0.075, down_freq=2.0, down_phase=0.0, gate_x=1.0,thick=0.2,offset=0.0,fade_rad=0.5,fade_z=0):
        """undulate the mesh with sin waves using the displacement values
        """
        down_freq=int(down_freq*gate_x)/gate_x
        up_freq=down_freq
        down_phase=0
        up_phase=0



        for vertex in self.mesh.vertices():
            vrt_attrs = self.mesh.vertex_attributes(vertex)
            orig_x = vrt_attrs["x"]
            orig_z = vrt_attrs["z"]
            x_new = vrt_attrs["x"] - vrt_attrs["x_disp"] #*0.99 # calc new x location
            z_new = vrt_attrs["z"] - vrt_attrs["z_disp"] #*0.99 # calc new x location

            fader=1.0
            dist_fade=math.sqrt((gate_x-x_new)*(gate_x-x_new)+(fade_z-z_new)*(fade_z-z_new))
            if dist_fade<fade_rad:
                fader=0.5+0.5*(1.0-math.cos((math.pi*dist_fade/fade_rad)))/2.0

            # z coordinate of the point normalized by the height of the wall
            rel_z = z_new/(self.z_size*self.elsize)

            # calculate y value of the up and down part of the wall at this x
            help1=1+thick*self.sin_wave(1, down_freq, down_phase, x_new)
            help2=1-thick*self.sin_wave(1, down_freq, down_phase, x_new)

            up_y = help1* self.sin_wave(up_amp, up_freq, up_phase, x_new)
            down_y = help1 * self.sin_wave(down_amp, down_freq, down_phase, x_new)

            # specify a y value as a linear blending of the top and the botom of the wall
            y_val = rel_z*up_y + (1.0-rel_z)*down_y

            # shift the y value in relation to the max amplitude:
            if len(list(self.mesh_offset.vertices())):
                #layer_amp = rel_z*up_amp + (1.0-rel_z)*down_amp
                #shift_val = down_amp - layer_amp
                #y_val = y_val + shift_val
                #self.mesh_offset.vertex_attribute(vertex, name='y', value=-y_val+2*down_amp)
                up_y_2 = help2* self.sin_wave(up_amp, up_freq, up_phase, x_new)
                down_y_2 = help2 * self.sin_wave(down_amp, down_freq, down_phase, x_new)
                y_val_2 = rel_z*up_y_2 + (1.0-rel_z)*down_y_2
                self.mesh_offset.vertex_attribute(vertex, name='y', value=fader*y_val_2-offset)

            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=fader*y_val)
                

    def create_mesh_displaced(self):
        """displace vertices with prior computed displacement vec from the gate points
        (this method is used just for preview)
        """

        self.mesh_displaced = self.mesh.copy()

        for vertex in self.mesh_displaced.vertices():
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')

            x_disp = self.mesh.vertex_attribute(vertex, name='x_disp')
            z_disp = self.mesh.vertex_attribute(vertex, name='z_disp')

            x_new = x_val - x_disp*0.99009 #for preview only, otherwise mesh cannot be drawn
            z_new = z_val - z_disp*0.99009 

            self.mesh_displaced.vertex_attribute(vertex, name='x', value=x_new)
            self.mesh_displaced.vertex_attribute(vertex, name='z', value=z_new)

    def make_gate(self, gate_size, gate_x, gate_type):
        """calculate the gate points, and their relative position to the gate center on floor
        """        
        for vertex in self.mesh.vertices():
            # get x and z coordinates from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')

            # calculate distance
            gate_rel_x = gate_x-x_val
            gate_rel_z = z_val

            #in order to avoid division by zero. should be cleaned up ASAP
            if gate_rel_x==0:
                gate_rel_x=0.000000001

            #this boolean will be set to True if vertex lies inside the gate 
            in_gate_vertex = False

            #for each gate tpye, the a corresponding function is called to check
            # to check if vertex lies inside the gate
            if gate_type=="circular":
                if self.in_gate_circular(gate_size, gate_rel_x, gate_rel_z):
                    in_gate_vertex = True
            if gate_type=="persian":
                if self.in_gate_persian(gate_size, gate_rel_x, gate_rel_z):
                    in_gate_vertex = True

            if in_gate_vertex:
                # store new vertex attribute with the distance to the gate center
                self.mesh.vertex_attribute(vertex,"x_disp", value=-gate_rel_x)
                self.mesh.vertex_attribute(vertex,"z_disp", value=gate_rel_z)
                #in_gate flag indicates that this vertex lies in the gate
                self.mesh.vertex_attribute(vertex,"in_gate", value=True)
                self.gate_points.append(vertex)

    def in_gate_circular(self, gate_size, gate_rel_x, gate_rel_z):
        """returns True if the relative coordintes (gate_rel_x,gate_rel_z) lie in
        a circular gate of size gate_size. otherwise false
        """   

        # calculate region (distance to gate center)
        dist = math.sqrt(gate_rel_x*gate_rel_x+gate_rel_z*gate_rel_z)
        if dist < gate_size:
            return True
        else:
            return False    
    
    def in_gate_triangular(self, gate_size, gate_rel_x, gate_rel_z):
        """returns True if the relative coordintes (gate_rel_x,gate_rel_z) lie in
        a triangular gate of size gate_size. otherwise false
        """   

        # calculate region (distance to gate center)
        if (gate_rel_z + abs(gate_rel_x)) < gate_size:
            return True
        else:
            return False    

    def in_gate_persian(self, gate_size, gate_rel_x, gate_rel_z):
        """returns Tre if the relative coordintes (gate_rel_x,gate_rel_z) lie in
        a "persian" gate of size gate_size. otherwise false
        """   

        # calculate regions
        # region 1 is the lower part (rectangular)
        region1 = gate_rel_z < gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0)
        # equation 1 is true if the x distance to gate center is smaller than gate width / 2
        equation1 = abs(gate_rel_x)-gate_size < self.elsize

        #region 2 is the midd section with smaller radii
        region2 = gate_rel_z > gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0) and gate_rel_z < gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2))
        equation2 = abs((gate_rel_x)**2+(gate_rel_z-gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0))**2)-gate_size*gate_size < self.elsize*self.elsize
        
        #region 3 is the top section with larger radii and the pointy peak
        region3 = gate_rel_z > gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)) #and gate_rel_z < gate_size*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(7/2))
        equation3 = abs((gate_rel_x+gate_rel_x*gate_size*math.sqrt(2)/(2*abs(gate_rel_x)))**2+(gate_rel_z-gate_size*(math.sqrt(2.5-math.sqrt(2.0))))**2)-4*gate_size*gate_size < self.elsize*self.elsize

        #checking if the point (coordinate) in question lies in  any of the regions
        if (region1 and equation1) or (region2 and equation2) or (region3 and equation3):
            return True
        else:
            return False

    def numpy_test_function(self, value):
        from compas.rpc import Proxy
        np = Proxy('numpy')
        linalg = Proxy('numpy.linalg')

        print("running numpy test function")

        a = np.array([[1, 2], [3, 5]])
        b = np.array([1, 2*value])
        x = linalg.solve(a, b)

        return x