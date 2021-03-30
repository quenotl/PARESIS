import bpy

from random import random
import numpy as np
from bpy import context
import json

membraneSize = 70.0 #(en dizaine de um)
sideSize = 150.0 #(en dizaine de um)
numberOfSpheres = 4000

# ---- Bottom Membrane  "ground" ----
bpy.ops.mesh.primitive_plane_add(size=membraneSize, align='WORLD', location=(0, 0, 0.1))
QSD = bpy.data.materials.new(name="QSD")
bpy.context.active_object.data.materials.append(QSD) #add the material to the object
bpy.context.object.active_material.diffuse_color = (0, 0, 1, 1) #Blue color

# ---- Adding physics to the membrane "ground" ----
bpy.ops.rigidbody.object_add()
bpy.context.object.rigid_body.type = 'PASSIVE'
bpy.context.object.rigid_body.collision_shape = 'MESH'
bpy.context.object.rigid_body.collision_margin = 0
bpy.context.object.rigid_body.friction = 0.046875
bpy.context.object.rigid_body.restitution = 0.046875



# ---- Lateral "wall" 1 ----
bpy.ops.mesh.primitive_plane_add(size=sideSize, enter_editmode=False, align='WORLD', location=(0, 0, 0))
bpy.context.object.location[0] = membraneSize/2.0
bpy.context.object.rotation_euler[1] = 1.5708
bpy.context.object.location[2] = sideSize/2
bpy.context.object.scale[1] = membraneSize/sideSize


# ---- Adding physics to the lateral "wall" 1 ----
bpy.ops.rigidbody.object_add()
bpy.context.object.rigid_body.type = 'PASSIVE'
bpy.context.object.rigid_body.collision_shape = 'MESH'
bpy.context.object.rigid_body.collision_margin = 0
bpy.context.object.rigid_body.friction = 0.046875
bpy.context.object.rigid_body.restitution = 0.046875
Side = bpy.data.materials.new(name="Side")
bpy.context.active_object.data.materials.append(Side) #add the material to the object
bpy.context.object.active_material.diffuse_color = (0, 0, 0, 0) #alpha = 0.0




# ---- Lateral "wall" 1 ----
bpy.ops.mesh.primitive_plane_add(size=sideSize, enter_editmode=False, align='WORLD', location=(0, 0, 0))
bpy.context.object.location[0] = -membraneSize/2.0
bpy.context.object.rotation_euler[1] = 1.5708
bpy.context.object.location[2] = sideSize/2
bpy.context.object.scale[1] = membraneSize/sideSize

# ---- Adding physics to the lateral "wall" 1 ----
bpy.ops.rigidbody.object_add()
bpy.context.object.rigid_body.type = 'PASSIVE'
bpy.context.object.rigid_body.collision_shape = 'MESH'
bpy.context.object.rigid_body.collision_margin = 0
bpy.context.object.rigid_body.friction = 0.046875
bpy.context.object.rigid_body.restitution = 0.046875
bpy.context.active_object.data.materials.append(Side) #add the material to the object
bpy.context.object.active_material.diffuse_color = (0, 0, 0, 0) #alpha = 0.0



# ---- Lateral "wall" 1 ----
bpy.ops.mesh.primitive_plane_add(size=sideSize, enter_editmode=False, align='WORLD', location=(0, 0, 0))
bpy.context.object.location[1] = membraneSize/2.0
bpy.context.object.rotation_euler[0] = 1.5708
bpy.context.object.location[2] = sideSize/2
bpy.context.object.scale[0] = membraneSize/sideSize

# ---- Adding physics to the lateral "wall" 1 ----
bpy.ops.rigidbody.object_add()
bpy.context.object.rigid_body.type = 'PASSIVE'
bpy.context.object.rigid_body.collision_shape = 'MESH'
bpy.context.object.rigid_body.collision_margin = 0
bpy.context.object.rigid_body.friction = 0.046875
bpy.context.object.rigid_body.restitution = 0.046875
bpy.context.active_object.data.materials.append(Side) #add the material to the object
bpy.context.object.active_material.diffuse_color = (0, 0, 0, 0) #alpha = 0.0



# ---- Lateral "wall" 1 ----
bpy.ops.mesh.primitive_plane_add(size=sideSize, enter_editmode=False, align='WORLD', location=(0, 0, 0))
bpy.context.object.location[1] = -membraneSize/2.0
bpy.context.object.rotation_euler[0] = 1.5708
bpy.context.object.location[2] = sideSize/2
bpy.context.object.scale[0] = membraneSize/sideSize

# ---- Adding physics to the lateral "wall" 1 ----
bpy.ops.rigidbody.object_add()
bpy.context.object.rigid_body.type = 'PASSIVE'
bpy.context.object.rigid_body.collision_shape = 'MESH'
bpy.context.object.rigid_body.collision_margin = 0
bpy.context.object.rigid_body.friction = 0.046875
bpy.context.object.rigid_body.restitution = 0.046875
bpy.context.active_object.data.materials.append(Side) #add the material to the object
bpy.context.object.active_material.diffuse_color = (0, 0, 0, 0) #alpha = 0.0


#Generating the spheres
for i in range(0, numberOfSpheres):
    
    #Progress bar (cuz it's fucking long)
    percent = float(i) * 100 / numberOfSpheres
    arrow   = '-' * int(percent/100 * 20 - 1) + '>'
    spaces  = ' ' * (20 - len(arrow))

    print('Simulation progress: [%s%s] %d %%' % (arrow, spaces, percent), end='\r')
    
    #Radius distribution (loi normale)
    mu = 6.1/10
    sigma = 2.28/10
    rayon = sigma * np.random.randn() + mu
    
    #Sphere's initial position (loi uniforme)
    minX = -membraneSize/2 * 0.9
    maxX = +membraneSize/2 * 0.9
    minY = -membraneSize/2 * 0.9
    maxY = +membraneSize/2 * 0.9

    x = minX+(maxX-minX)*random()
    y = minY+(maxY-minY)*random()

    z = 130 - (100*rayon) + 10*random() #Z value depends on the radius, the smaller, the higher

    if z < 0: #Just in case a sphere appears under the ground
        z = rayon + 0.1
    
    

    if rayon > 0: #Only taking positive radia into account
        bpy.ops.mesh.primitive_uv_sphere_add(radius = rayon, location=(x,y,z))
        bpy.ops.rigidbody.object_add()

        #Adding physics to the spheres
        bpy.context.object.rigid_body.mass = 4/3*3.14*rayon*rayon*rayon #mass depends on the size
        bpy.context.object.rigid_body.collision_shape = 'SPHERE'
        bpy.context.object.rigid_body.friction = 1
        bpy.context.object.rigid_body.collision_margin = 0
        bpy.context.object.rigid_body.use_margin = True
        bpy.context.object.rigid_body.linear_damping = 0.35
        bpy.context.object.rigid_body.angular_damping = 0.6

        #Attibuting a color depending on the size (red = BIG, green = SMALL)
        distColor = rayon/1.2
        Color = bpy.data.materials.new(name="Color" + str(i))
        bpy.context.active_object.data.materials.append(Color) #add the material to the object
        bpy.context.object.active_material.diffuse_color = (distColor, 1-distColor, 0, 1)

#Runnin simulation until frame 229
for i in range(0,265):
    bpy.context.scene.frame_set(i)
    bpy.context.view_layer.update()

scene = context.scene
spheres = [o for o in scene.objects if o.name.startswith("Sphere")]

#Checking sphere's final positions
sphereList = []

for s in spheres:
    o = s.matrix_world.to_translation()
    radius = sum(s.dimensions) / 6 # sphere radius estimation
    
    difference = radius - s.matrix_world.translation[2]
    if difference < 0.95*radius and difference > -0.95*radius:
        sphereList.append([s.matrix_world.translation[0]*10, s.matrix_world.translation[1]*10, radius*10])
    else:
        Transparent = bpy.data.materials.new(name="Transparent" )
        s.data.materials.append(Transparent) #add the material to the object
        s.active_material.diffuse_color = (0, 0, 0, 0)
      
        

#Saving results
with open("/Users/quenot/Documents/Simulations/rdmSphere_surf700x700um_meanRad6p1um_sigma2p28_b.txt", "a") as fp:
    json.dump(sphereList, fp)