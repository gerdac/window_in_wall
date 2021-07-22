from compas.datastructures import Mesh

my_mesh = Mesh()

u = 10
v = 5

for i in range(u):
    for j in range(v):
        my_mesh.add_vertex(x=i, y=j, z=0)


for i in range(u-1):
    for j in range(v-1):
        a = i * v + j
        b = a + v
        c = b + 1
        d = a + 1
        my_mesh.add_face([a,b,c,d])


