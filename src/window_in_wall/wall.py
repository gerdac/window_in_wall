import math
import compas

from compas.datastructures import Mesh


class Wall(object):
    """ Wall generator

     """

    def __init__(self):
        """ constructor

        """
        self.mesh = Mesh()
        self.gate_points = []
        self.printTopol = Mesh ()
        self.printFrames = []


    def create_quad_mesh_from_dimensions(self, width=3.5, height=3.2, elsize=0.025):
        """create a regular quad mesh from height and width and edge length in the xz-plane
        """
        
        self.elsize = elsize

        self.x_size = int(width / elsize) + 1
        self.z_size = int(height / elsize) + 1

        # create vertices
        for i in range(self.x_size):
            for j in range(self.z_size):
                x = i * self.elsize
                y = 0
                z = j * self.elsize
                glob_id = i*self.z_size+j
                self.mesh.add_vertex(x = x, y = 0, z = z, i=i, j=j, glob_id=glob_id, x_disp=0, z_disp=0)

        # create faces
        for i in range(self.x_size - 1):
            for j in range(self.z_size - 1):
                a = i * self.z_size + j
                b = a + self.z_size
                c = b + 1
                d = a + 1
                self.mesh.add_face([a, b, c, d])
    
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
    
    def undulate_with_displ(self, up_amp=0.075, up_freq=2.0, up_phase=0.0,down_amp=0.075, down_freq=2.0, down_phase=0.0):
        """undulate the mesh with sin waves using the displacement values
        """
        for vertex in self.mesh.vertices():
            
            x_val = self.mesh.vertex_attribute(vertex, name='x') # get x coordinate from vertex
            z_val = self.mesh.vertex_attribute(vertex, name='z')
            x_disp = self.mesh.vertex_attribute(vertex, name='x_disp') # get displacement value fomr vertex
            x_new = x_val - x_disp #*0.99 # calc new x location
            z_disp = self.mesh.vertex_attribute(vertex, name='z_disp') # get displacement value fomr vertex
            z_new = z_val - z_disp #*0.99 # calc new x location

            # z coordinate of the point normalized by the height of the wall
            rel_z=z_new/(self.z_size*self.elsize)
            # calculate y value of the up and down part of the wall at this x
            up_y = self.sin_wave(up_amp, up_freq, up_phase, x_new)
            down_y = self.sin_wave(down_amp, down_freq, down_phase, x_new)
            # specify a y value as a linear blending of the top and the botom of the wall
            y_val = rel_z*up_y + (1.0-rel_z)*down_y

            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=y_val)

    def displace_vertices(self):
        """displace vertices with prior computed displacement vec from the gate points
        (this method is used just for preview)
        """

        for vertex in self.mesh.vertices():
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')

            x_disp = self.mesh.vertex_attribute(vertex, name='x_disp')
            z_disp = self.mesh.vertex_attribute(vertex, name='z_disp')

            x_new = x_val - x_disp*0.99009 #for preview only, otherwise mesh cannot be drawn
            z_new = z_val - z_disp*0.99009 

            self.mesh.vertex_attribute(vertex, name='x', value=x_new)
            self.mesh.vertex_attribute(vertex, name='z', value=z_new)

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
                if self.in_gate_circular(gate_size, gate_rel_x,gate_rel_z):
                    in_gate_vertex = True
            if gate_type=="persian":
                if self.in_gate_persian(gate_size, gate_rel_x,gate_rel_z):
                    in_gate_vertex = True

            if in_gate_vertex:
                # store new vertex attribute with the distance to the gate center
                self.mesh.vertex_attribute(vertex,"x_disp", value=-gate_rel_x)
                self.mesh.vertex_attribute(vertex,"z_disp", value=gate_rel_z)
                #in_gate flag indicates that this vertex lies in the gate
                self.mesh.vertex_attribute(vertex,"in_gate", value=True)
                self.gate_points.append(vertex)

    def in_gate_circular(self,gate_size, gate_rel_x,gate_rel_z):
        """returns True if the relative coordintes (gate_rel_x,gate_rel_z) lie in
        a circular gate of size gate_size. otherwise false
        """   

        # calculate region (distance to gate center)
        dist = math.sqrt(gate_rel_x*gate_rel_x+gate_rel_z*gate_rel_z)
        if dist < gate_size:
            return True
        else:
            return False    
    
    def in_gate_persian(self,gate_size, gate_rel_x,gate_rel_z):
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