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
                self.mesh.add_vertex(x = x, y = 0, z = z, i=i, j=j)

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

    def undulate(self, amp=0.075, freq=2.0, phase=0.0):
        """undulate the mesh with sin waves
        """
        for vertex in self.mesh.vertices():
            # get x coordinate from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            # calculate y value offset from x
            y_val = self.sin_wave(amp, freq, phase, x_val)
            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=y_val)
    
    def gate_circular(self, gate_radius=1.0, gate_x=1.75):
        """calculate the gate points of a circular gate
        """
        gate_radius = gate_radius * 1.0

        for vertex in self.mesh.vertices():
            # get x and z coordinates from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')

            # calculate distance
            xx = gate_x-x_val
            zz = z_val

            # calculate distance to the x location of the gate
            dist = math.sqrt((x_val-gate_x)*(x_val-gate_x)+z_val*z_val)
            if abs(dist-gate_radius)<self.elsize:
                # store new vertex attribute with the distance to the gate center
                self.mesh.vertex_attribute(vertex, name='gate_point_value', value=[-xx,zz])
                i = self.mesh.vertex_attribute(vertex, "i")
                j = self.mesh.vertex_attribute(vertex, "j")
                self.gate_points.append([i,j,-xx,zz])
    
    def gate_persian(self, gate_radius=1.0, gate_x=1.75):
        """calculate the gate points of a persian gate
        """
        gate_radius = gate_radius * 1.0

        for vertex in self.mesh.vertices():
            # get x and z coordinates from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            z_val = self.mesh.vertex_attribute(vertex, name='z')

            # calcluate distance
            xx = gate_x-x_val
            if xx==0:
                xx=0.0000001
            zz = z_val

            # calculate regions
            region1 = zz < gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0)
            equation1 = abs(xx)-gate_radius < self.elsize

            region2 = zz > gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0) and zz < gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2))
            equation2 = abs((xx)**2+(zz-gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)/2.0))**2)-gate_radius*gate_radius < self.elsize*self.elsize
            
            region3 = zz > gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(2)) and zz < gate_radius*(math.sqrt(2.5-math.sqrt(2.0))+math.sqrt(7/2))
            equation3 = abs((xx+xx*gate_radius*math.sqrt(2)/(2*abs(xx)))**2+(zz-gate_radius*(math.sqrt(2.5-math.sqrt(2.0))))**2)-4*gate_radius*gate_radius < self.elsize*self.elsize

            if (region1 and equation1) or (region2 and equation2) or (region3 and equation3):
                # store new vertex attribute with the distance to the gate center
                self.mesh.vertex_attribute(vertex, name='gate_point_value', value=[-xx,zz])
                i = self.mesh.vertex_attribute(vertex, "i")
                j = self.mesh.vertex_attribute(vertex, "j")
                self.gate_points.append([i,j,-xx,zz])

    def numpy_test_function(self, value):
        from compas.rpc import Proxy
        np = Proxy('numpy')
        linalg = Proxy('numpy.linalg')

        print("running numpy test function")

        a = np.array([[1, 2], [3, 5]])
        b = np.array([1, 2*value])
        x = linalg.solve(a, b)

        return x