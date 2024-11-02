import numpy as np
from matplotlib import pyplot as plt
from src.panel_method import Airfoil, SurfacePanelMesh, WakePanelMesh, AeroPanelMesh, PanelMethod



def simulate_steady_flow_around_an_Airfoil(
    airfoil_name="naca0012 sharp",
    chord_length=5,
    leading_edge_location=(3, 1),
    velocity=1,
    angle_of_attack=0,
    num_airfoil_panels=30,
    wake_length_in_chords=15,
    num_wake_panels=1,
    kutta_vertex_id=0,
    wake_relaxation_iters=0
):
       
    airfoil = Airfoil(
        name=airfoil_name,
        chord_length=chord_length,
        num_points=num_airfoil_panels+1
    )
    vertex = airfoil.coords[0:-1]   
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    
    surface_mesh = SurfacePanelMesh(
        vertex=vertex,
        face=face,
        CCW=True,
        kutta_vertex_id=kutta_vertex_id
    )
    
    wake_mesh = WakePanelMesh(
        surface_mesh=surface_mesh,
        length = wake_length_in_chords * airfoil.chord,
        num_faces= num_wake_panels,
        surface_fixed=True
    )
    
    panel_method = PanelMethod(
        mesh=AeroPanelMesh(
            surface_mesh=surface_mesh,
            wake_mesh=wake_mesh
        )
    )
    
    
    panel_method.set_BodyFixedFrame_origin(
        leading_edge_location[0], leading_edge_location[1]
    )
    
    panel_method.set_V_fs(angle_of_attack, velocity)
    
    
    panel_method.set_BodyFixedFrame_origin(
        leading_edge_location[0], leading_edge_location[1]
    )
    
    panel_method.set_V_fs(angle_of_attack, velocity)
    # panel_method.set_V_inf(angle_of_attack, velocity)
    
    
    panel_method.mesh.display(
        BodyFixed_FrameOfReference=True,
        display_normals=True
    )

    panel_method.mesh.display(
        BodyFixed_FrameOfReference=False,
        display_normals=True
    )
    
    
    panel_method.solve(steady_state=True, iters=wake_relaxation_iters)
    
    print("CF = ", panel_method.AerodynamicForce(CharLength=airfoil.chord))
    print("CL = ", panel_method.LiftForce(CharLength=airfoil.chord).norm())
    print(
        "CD = ", panel_method.InducedDragForce(CharLength=airfoil.chord).norm()
    )
    
    print("r_cop = ", panel_method.Center_of_Pressure(CharLength=airfoil.chord))
    
    panel_method.mesh.display(
        BodyFixed_FrameOfReference=False, display_normals=False
    )
    panel_method.display_forces()
    
    plot_pressure_coefficient(panel_method, airfoil)
    plot_fields(panel_method, airfoil)   

def simulate_unsteady_flow_around_an_Airfoil(
    airfoil_name="naca0012 sharp",
    chord_length=5,
    leading_edge_location=(15, 2),
    velocity=1,
    angle_of_attack=10,
    num_airfoil_panels=30,
    kutta_vertex_id=0,
    num_time_steps=40
):
    
    airfoil = Airfoil(
        name=airfoil_name,
        chord_length=chord_length,
        num_points=num_airfoil_panels+1
    )
    vertex = airfoil.coords[0:-1]   
    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])
    
    
    surface_mesh = SurfacePanelMesh(
        vertex=vertex,
        face=face,
        CCW=True,
        kutta_vertex_id=kutta_vertex_id
    )
    
    surface_mesh.set_BodyFixedFrame_origin(
        leading_edge_location[0], leading_edge_location[1]
    )
    
    surface_mesh.set_BodyFixedFrame_orientation(-angle_of_attack)
    
    wake_mesh = WakePanelMesh(
        surface_mesh=surface_mesh,
        length = 0,
        num_faces= 0,
        surface_fixed=False
    )
    
    panel_method = PanelMethod(
        mesh=AeroPanelMesh(
            surface_mesh=surface_mesh,
            wake_mesh=wake_mesh
        )
    )
    
            
    panel_method.set_BodyFixedFrame_origin_velocity(Vo_x=-velocity, Vo_y=0)
    
    panel_method.mesh.display(
        BodyFixed_FrameOfReference=False, display_normals=True
    )
        
    panel_method.solve(steady_state=False, iters=num_time_steps)
    
    print("CF = ", panel_method.AerodynamicForce(CharLength=airfoil.chord))
    print("CL = ", panel_method.LiftForce(CharLength=airfoil.chord).norm())
    print(
        "CD = ", panel_method.InducedDragForce(CharLength=airfoil.chord).norm()
    )
    
    print("r_cop = ", panel_method.Center_of_Pressure(CharLength=airfoil.chord))
    
    panel_method.mesh.display(
        BodyFixed_FrameOfReference=False, display_normals=False
    )
    panel_method.display_forces()
    
    plot_pressure_coefficient(panel_method, airfoil)
    
    plot_fields(panel_method, airfoil)
  
def meshgrid(panel_method:PanelMethod, length:float, num:int):

    # only works good velocity field and contour plot
    # it doesnt work with streamplot
    
    num_s = len(panel_method.surface.panel)
    num_wake = len(panel_method.wake.panel)
    
    X=np.zeros((num_s + 2 * num_wake, num))
    Y = np.zeros_like(X)
    
    
    for j in range(num):
          
        for surface_panel in panel_method.surface.panel:
            i = surface_panel.id
            r = surface_panel.r_cp + (j+1)*length * surface_panel.e_n
            r = panel_method.surface.ro + r.transform(panel_method.surface.A.T)
            X[i][j], Y[i][j] = r.x, r.y
            pass
        
        for wake_panel in panel_method.wake.panel:
            
            i_up = wake_panel.id + num_s
            r_up = wake_panel.r_cp + (j+1)*length * wake_panel.e_n
            
            i_down = wake_panel.id + num_s + num_wake
            r_down = wake_panel.r_cp - (j+1)*length * wake_panel.e_n
            
            r_up = panel_method.wake.ro + r_up.transform(panel_method.wake.A.T)
            r_down = panel_method.wake.ro + r_down.transform(panel_method.wake.A.T)
            
            if panel_method.wake.surface_fixed:
                                
                r_up = (
                    panel_method.surface.ro
                    + r_up.transform(panel_method.surface.A.T)
                )
                
                r_down = (
                    panel_method.surface.ro
                    + r_down.transform(panel_method.surface.A.T)
                )      

            X[i_up][j], Y[i_up][j] = r_up.x, r_up.y
            X[i_down][j], Y[i_down][j] = r_down.x, r_down.y
            
    
    return X, Y
    
def plot_pressure_coefficient(panel_method, airfoil):
    plt.plot(
        [panel.r_cp.x/airfoil.chord for panel in panel_method.surface.panel],
        [panel.Cp for panel in panel_method.surface.panel],
        'ks--', markerfacecolor='r', label='Panel Method'
    )
    plt.xlabel("x/c")
    plt.ylabel("Cp")
    plt.title("Chordwise Pressure Coefficient Distribution")
    plt.legend()
    plt.grid()
    plt.gca().invert_yaxis()
    plt.show()
    
def plot_fields(panel_method, airfoil):
    
    X, Y = np.meshgrid(
        np.linspace(
            start=panel_method.surface.ro.x - airfoil.chord,
            stop=panel_method.surface.ro.x + 2 * airfoil.chord,
            num=100
        ),
        np.linspace(
            start=panel_method.surface.ro.y - airfoil.chord/2,
            stop=panel_method.surface.ro.y + airfoil.chord/2,
            num=50
        ),
        indexing='ij'
    )
        
    
    panel_method.display_contour(X, Y, BodyFixed_FrameOfReference=False)
    panel_method.display_velocity_field(X, Y, BodyFixed_FrameOfReference=False)
    panel_method.display_streamlines(X, Y, BodyFixed_FrameOfReference=False)


def test_PanelMethod_SteadyState_rigidWake():
    
    simulate_steady_flow_around_an_Airfoil(
        airfoil_name="naca0012 sharp",
        chord_length=5,
        leading_edge_location=(3, 1),
        velocity=1,
        angle_of_attack=10,
        num_airfoil_panels=30,
        wake_length_in_chords=15,
        num_wake_panels=20,
        kutta_vertex_id=0,
        wake_relaxation_iters=0
    )
      
def test_PanelMethod_SteadyState_iterativeWake():
    
    simulate_steady_flow_around_an_Airfoil(
        airfoil_name="naca0012 sharp",
        chord_length=5,
        leading_edge_location=(3, 1),
        velocity=1,
        angle_of_attack=10,
        num_airfoil_panels=30,
        wake_length_in_chords=15,
        num_wake_panels=20,
        kutta_vertex_id=0,
        wake_relaxation_iters=10
    )    

def test_PanelMethod_Unsteady():
    simulate_unsteady_flow_around_an_Airfoil(
        airfoil_name="naca0012 sharp",
        chord_length=5,
        leading_edge_location=(15, 2),
        velocity=1,
        angle_of_attack=10,
        num_airfoil_panels=30,
        kutta_vertex_id=0,
        num_time_steps=40
    )
    