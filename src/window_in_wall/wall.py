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


    def create_quad_mesh_from_dimensions(self, width, height, elsize):
        """create a regular quad mesh from height and width and edge length in the xz-plane
        """

        x_size = int(width / elsize) + 1
        z_size = int(height / elsize) + 1

        # create vertices
        for i in range(x_size):
            for j in range(z_size):
                x = i * elsize
                y = 0
                z = j * elsize
                self.mesh.add_vertex(x = x, y = 0, z = z)

        # create faces
        for i in range(x_size - 1):
            for j in range(z_size - 1):
                a = i * z_size + j
                b = a + z_size
                c = b + 1
                d = a + 1
                self.mesh.add_face([a, b, c, d])
    
    def sin_wave(self, amp, freq, phase, value):
        return(amp * math.sin(2.0 * math.pi * freq * value + phase))

    def undulate(self, amp, freq, phase):
        """undulate the mesh with sin waves
        """
        for vertex in self.mesh.vertices():
            # get x coordinate from vertex
            x_val = self.mesh.vertex_attribute(vertex, name='x')
            # calculate y value offset from x
            y_val = self.sin_wave(amp, freq, phase, x_val)
            # set y coordinate from vertex
            self.mesh.vertex_attribute(vertex, name='y', value=y_val)