import numpy as  np
from src.myMath import Vector
from src.panel_method import Circle, Airfoil, Mesh, SurfaceMesh, PanelMesh, AeroPanelMesh, WakeMesh, WakePanelMesh, SurfacePanelMesh, AeroMesh


def test_mesh():
    
    circle = Circle(name="circle", center=(0, 0), radius=1, num_points=11)
    vertex = circle.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    mesh = Mesh(vertex=vertex, face=face)
    
    mesh.set_BodyFixedFrame_origin(xo=2*circle.radius, yo=2*circle.radius)
    mesh.set_BodyFixedFrame_orientation(theta_z=45)
    mesh.set_BodyFixedFrame_origin_velocity(Vo_x=circle.radius, Vo_y=0)
    mesh.set_BodyFixedFrame_angular_velocity(omega_z=-45)

    mesh.display(BodyFixed_FrameOfReference=False)    
    mesh.display(BodyFixed_FrameOfReference=True)
    
    mesh.move_BodyFixedFrame(dt=1)
    
    mesh.display(BodyFixed_FrameOfReference=False)    
    mesh.display(BodyFixed_FrameOfReference=True)
    
    mesh.move_BodyFixedFrame(dt=1)
    
    mesh.display(BodyFixed_FrameOfReference=False)    
    mesh.display(BodyFixed_FrameOfReference=True)
    
def test_surface_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=2, num_points=5)
    vertex = airfoil.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_mesh = SurfaceMesh(vertex=vertex, face=face, kutta_vertex_id=None)
    surface_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord/2, yo=airfoil.chord/2
    )
    surface_mesh.set_BodyFixedFrame_orientation(theta_z=-20)
    
    surface_mesh.display(BodyFixed_FrameOfReference=False)
    
    print("kutta vertex id: " + str(surface_mesh.kutta_vertex_id))
    print("bottom kuta face id: " + str(surface_mesh.bottom_kutta_face_id))
    print("top kutta face id: " + str(surface_mesh.top_kutta_face_id))
    print("adjacency matrix =")
    print(surface_mesh.adjacency_matrix)
    
    id=0
    print(
        "\nface with id: " + str(id) 
        + " has adjacent faces with id's: " + str(surface_mesh.give_adjacent_faces(face_id=id))
    )
    
    surface_mesh = SurfaceMesh(vertex=vertex, face=face, kutta_vertex_id=0)
    
    print("\nkutta vertex id: " + str(surface_mesh.kutta_vertex_id))
    print("bottom kuta face id: " + str(surface_mesh.bottom_kutta_face_id))
    print("top kutta face id: " + str(surface_mesh.top_kutta_face_id))
    print("adjacency matrix =")
    print(surface_mesh.adjacency_matrix)
    
    id=surface_mesh.top_kutta_face_id
    print(
        "\nface with id: " + str(id) 
        + " has an adjacent face with id: " + str(surface_mesh.give_adjacent_faces(face_id=id))
    )
    
    id=surface_mesh.bottom_kutta_face_id
    print(
        "\nface with id: " + str(id) 
        + " has an adjacent face with id: " + str(surface_mesh.give_adjacent_faces(face_id=id))
    )

def test_wake_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=2, num_points=5)
    vertex = airfoil.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_mesh = SurfaceMesh(vertex=vertex, face=face, kutta_vertex_id=0)
    surface_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord, yo=airfoil.chord
    )
    
    
    wake_mesh = WakeMesh(
        surface_mesh=surface_mesh,
        length= 5 * airfoil.chord,
        num_faces=1,
        surface_fixed=True
    )
    wake_mesh.set_BodyFixedFrame_orientation(15)
    wake_mesh.display(BodyFixed_FrameOfReference=False)
    wake_mesh.display(BodyFixed_FrameOfReference=True)
    
    wake_mesh = WakeMesh(
        surface_mesh=surface_mesh,
        length= 5 * airfoil.chord,
        num_faces=1,
        surface_fixed=False
    )
    wake_mesh.set_BodyFixedFrame_orientation(15)
    wake_mesh.display(BodyFixed_FrameOfReference=False)
    wake_mesh.display(BodyFixed_FrameOfReference=True)

def test_panel_mesh():
    
    circle = Circle(name="circle", center=(0, 0), radius=5, num_points=11)
    vertex = circle.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    
    panel_mesh = PanelMesh(vertex=vertex, face=face, CCW=True)
    
    panel_mesh.set_BodyFixedFrame_origin(xo=2*circle.radius, yo=2*circle.radius)
    panel_mesh.set_BodyFixedFrame_orientation(theta_z=45)
    panel_mesh.set_BodyFixedFrame_origin_velocity(Vo_x=circle.radius, Vo_y=0)
    panel_mesh.set_BodyFixedFrame_angular_velocity(omega_z=-45)

    panel_mesh.display(
        BodyFixed_FrameOfReference=False, display_normals=True
    )    
    panel_mesh.display(BodyFixed_FrameOfReference=True, display_normals=True
    )
    
    panel_mesh.move_BodyFixedFrame(dt=1)
    
    panel_mesh.display(BodyFixed_FrameOfReference=False, display_normals=True
    )    
    panel_mesh.display(BodyFixed_FrameOfReference=True, display_normals=True
    )
    
    panel_mesh.move_BodyFixedFrame(dt=1)
    
    panel_mesh.display(BodyFixed_FrameOfReference=False, display_normals=True
    )   
    panel_mesh.display(BodyFixed_FrameOfReference=True, display_normals=True
    )

def test_surface_panel_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=20, num_points=5)
    vertex = airfoil.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_panel_mesh = SurfacePanelMesh(
        vertex=vertex, face=face, CCW=True, kutta_vertex_id=None
    )
    surface_panel_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord/2, yo=airfoil.chord/2
    )
    surface_panel_mesh.set_BodyFixedFrame_orientation(theta_z=-20)
    
    surface_panel_mesh.display(
        BodyFixed_FrameOfReference=False, display_normals=True
    )
    
    print("kutta vertex id: " + str(surface_panel_mesh.kutta_vertex_id))
    print("bottom kuta face id: " + str(surface_panel_mesh.bottom_kutta_face_id))
    print("top kutta face id: " + str(surface_panel_mesh.top_kutta_face_id))
    print("adjacency matrix =")
    print(surface_panel_mesh.adjacency_matrix)
    
    id=0
    adjacent_panels = surface_panel_mesh.give_adjacent_panels(
        surface_panel_mesh.panel[id]
    )
    adjacent_panels_ids = [panel.id for panel in adjacent_panels]
    
    print(
        "\nface with id: " + str(id) 
        + " has adjacent panels with id's: " + str(adjacent_panels_ids)
    )
    
    surface_panel_mesh = SurfacePanelMesh(
        vertex=vertex, face=face, CCW=True, kutta_vertex_id=0
    )
    
    print("\nkutta vertex id: " + str(surface_panel_mesh.kutta_vertex_id))
    print("bottom kuta face id: " + str(surface_panel_mesh.bottom_kutta_face_id))
    print("top kutta face id: " + str(surface_panel_mesh.top_kutta_face_id))
    print("adjacency matrix =")
    print(surface_panel_mesh.adjacency_matrix)
    
    id=surface_panel_mesh.top_kutta_face_id
    adjacent_panels = surface_panel_mesh.give_adjacent_panels(
        surface_panel_mesh.panel[id]
    )
    adjacent_panels_ids = [panel.id for panel in adjacent_panels]
    
    print(
        "\npanel with id: " + str(id) 
        + " has an adjacent panel with id: " + str(adjacent_panels_ids)
    )
    
    id=surface_panel_mesh.bottom_kutta_face_id
    adjacent_panels = surface_panel_mesh.give_adjacent_panels(
        surface_panel_mesh.panel[id]
    )
    adjacent_panels_ids = [panel.id for panel in adjacent_panels]
    
    print(
        "\npanel with id: " + str(id) 
        + " has an adjacent panel with id: " + str(adjacent_panels_ids)
    )

def test_wake_panel_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=2, num_points=5)
    vertex = airfoil.coords[0:-1]
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_panel_mesh = SurfacePanelMesh(
        vertex=vertex, face=face, CCW=True, kutta_vertex_id=0
    )
    surface_panel_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord, yo=airfoil.chord
    )
    
    
    wake_panel_mesh = WakePanelMesh(
        surface_mesh=surface_panel_mesh,
        length = 10 * airfoil.chord,
        num_faces =5, 
        surface_fixed=True
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True
    )
    
    wake_panel_mesh.move_vertex(vertex_id=0, dr=Vector(0, airfoil.chord/5, 0))
    
    wake_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True
    )
    
    
    wake_panel_mesh = WakePanelMesh(
        surface_mesh=surface_panel_mesh,
        length = 10 * airfoil.chord,
        num_faces =5, 
        surface_fixed=False
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True
    )
    
    wake_panel_mesh.move_vertex(vertex_id=0, dr=Vector(0, airfoil.chord/5, 0))
    
    wake_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True
    )
    wake_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True
    )
    
def test_aero_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=2, num_points=20)

    vertex = airfoil.coords[0:-1]
    
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_mesh = SurfaceMesh(vertex=vertex, face=face, kutta_vertex_id=0)
    surface_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord/2, yo=airfoil.chord/2
    )
    surface_mesh.set_BodyFixedFrame_orientation(theta_z=-20)
    
    wake_mesh = WakeMesh(
        surface_mesh=surface_mesh,
        length= 3 * airfoil.chord,
        num_faces=1,
        surface_fixed=True
    )
    wake_mesh.set_BodyFixedFrame_orientation(20)
    
    aero_mesh = AeroMesh(surface_mesh, wake_mesh)
    
    aero_mesh.display(BodyFixed_FrameOfReference=False)
    
    aero_mesh.display(BodyFixed_FrameOfReference=True)
    
    wake_mesh = WakeMesh(
        surface_mesh=surface_mesh,
        length= 3 * airfoil.chord,
        num_faces=1,
        surface_fixed=False
    )
    wake_mesh.set_BodyFixedFrame_orientation(20)
    
    aero_mesh = AeroMesh(surface_mesh, wake_mesh)
    
    aero_mesh.display(BodyFixed_FrameOfReference=False)
    
    aero_mesh.display(BodyFixed_FrameOfReference=True)

def test_aero_panel_mesh():
    
    airfoil = Airfoil("naca0012 sharp", chord_length=20, num_points=5)

    vertex = airfoil.coords[0:-1]
    
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_panel_mesh = SurfacePanelMesh(
        vertex=vertex, face=face, CCW=True, kutta_vertex_id=0
    )
    surface_panel_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord/2, yo=airfoil.chord/2
    )
    
    surface_panel_mesh.set_BodyFixedFrame_orientation(-20)
    
    
    wake_panel_mesh = WakePanelMesh(
        surface_mesh=surface_panel_mesh,
        length = 2 * airfoil.chord,
        num_faces =3, 
        surface_fixed=True
    )
    
    wake_panel_mesh.set_BodyFixedFrame_orientation(20)
    
    aero_panel_mesh = AeroPanelMesh(surface_panel_mesh, wake_panel_mesh)
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    aero_panel_mesh.wake_mesh.move_vertex(
        vertex_id=len(aero_panel_mesh.wake_mesh.vertex)-1,
        dr=Vector(0, airfoil.chord/2, 0)
    )
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    
    wake_panel_mesh = WakePanelMesh(
        surface_mesh=surface_panel_mesh,
        length = 2 * airfoil.chord,
        num_faces = 3, 
        surface_fixed=False
    )
    wake_panel_mesh.set_BodyFixedFrame_orientation(20)
    
    aero_panel_mesh = AeroPanelMesh(surface_panel_mesh, wake_panel_mesh)
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    aero_panel_mesh.wake_mesh.move_vertex(
        vertex_id=0, dr=Vector(0, airfoil.chord/2, 0)
    )
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    pass  

def test_wake_kinematics():
    airfoil = Airfoil("naca0012 sharp", chord_length=20, num_points=5)

    vertex = airfoil.coords[0:-1]
    
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    surface_panel_mesh = SurfacePanelMesh(
        vertex=vertex, face=face, CCW=True, kutta_vertex_id=0
    )
    surface_panel_mesh.set_BodyFixedFrame_origin(
        xo=airfoil.chord/2, yo=airfoil.chord/2
    )
    
    surface_panel_mesh.set_BodyFixedFrame_angular_velocity(15)
    
    wake_panel_mesh = WakePanelMesh(
        surface_mesh=surface_panel_mesh,
        length = 1 * airfoil.chord,
        num_faces = 2, 
        surface_fixed=False
    )
    
    
    
    aero_panel_mesh = AeroPanelMesh(surface_panel_mesh, wake_panel_mesh)
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    for i in range(6):
        aero_panel_mesh.move_rigid_body(dt=1)
        aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    surface_panel_mesh.set_BodyFixedFrame_origin(
        xo=3*airfoil.chord, yo=airfoil.chord
    )
    surface_panel_mesh.set_BodyFixedFrame_orientation(-15)
    surface_panel_mesh.set_BodyFixedFrame_angular_velocity(0)
    surface_panel_mesh.set_BodyFixedFrame_origin_velocity(
        Vo_x=-airfoil.chord/2, Vo_y=0)
    wake_panel_mesh = WakePanelMesh(surface_panel_mesh, 0, 0, False)
    
    V_inf = Vector(airfoil.chord/2, 0, 0)
    wake_panel_mesh.set_BodyFixedFrame_origin_velocity(
        Vo_x=V_inf.x, Vo_y=V_inf.y
    )
      
    aero_panel_mesh = AeroPanelMesh(surface_panel_mesh, wake_panel_mesh)
    aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
    
    for i in range(5):
        aero_panel_mesh.move_rigid_body(dt=1)
        aero_panel_mesh.display(BodyFixed_FrameOfReference=False,
                            display_normals=True)
        aero_panel_mesh.wake_mesh.panel[-1].mu=1
    
    aero_panel_mesh = AeroPanelMesh(surface_panel_mesh, wake_panel_mesh)
    aero_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True)
    
    for i in range(5):
        aero_panel_mesh.move_rigid_body(dt=1)
        aero_panel_mesh.display(BodyFixed_FrameOfReference=True,
                            display_normals=True)
   