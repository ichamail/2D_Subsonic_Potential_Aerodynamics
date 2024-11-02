import numpy as np
from matplotlib import pyplot as plt
from src.panel_method import Circle, SurfacePanelMesh, BoundaryElementMethod

def simulate_flow_around_a_2D_circular_object(
    radius=5, center=(0, 0),velocity=1, angle_of_attack=0, num_panels=10
):

    circle = Circle(
        name="cirlce", center=(0, 0), radius=radius, num_points=num_panels+1
    )
    vertex = circle.coords[0:-1]

    face = np.array([[i, (i+1)%len(vertex)] for i in range(len(vertex))])

    surface_mesh = SurfacePanelMesh(
        vertex=vertex,
        face=face,
        CCW=True,
        kutta_vertex_id=None
    )

    panel_method = BoundaryElementMethod(surface_mesh)

    panel_method.set_BodyFixedFrame_origin(center[0], center[1])
    panel_method.set_V_fs(angle_of_attack=angle_of_attack, magnitude=velocity)
    # panel_method.set_V_inf(angle_of_attack=angle_of_attack, magnitude=velocity)

    panel_method.surface.display(
        BodyFixed_FrameOfReference=True,
        display_normals=True
    )


    panel_method.surface.display(
        BodyFixed_FrameOfReference=False,
        display_normals=True
    )

    
    panel_method.solve()
    
    plot_pressure_coefficient(panel_method, angle_of_attack)
    plot_fields(panel_method, circle)
    
def plot_pressure_coefficient(panel_method, angle_of_attack:float):
    
    analytical_theta = np.linspace(
        -np.pi + np.deg2rad(angle_of_attack),
        np.pi + np.deg2rad(angle_of_attack),
        200
    )
    analytical_cp = 1 - 4 * np.sin(analytical_theta)**2
    plt.plot(analytical_theta*(180/np.pi), analytical_cp ,'g-',
            label='Analytical')

    plt.plot(
        [
            np.arctan2(panel.r_cp.y, panel.r_cp.x) * (180/np.pi) 
            + angle_of_attack
            for panel in panel_method.surface.panel
        ],
        [panel.Cp for panel in panel_method.surface.panel],
        'ks', markerfacecolor='r', label='Panel Method'
    )

    plt.xlabel('Angle [deg]')
    plt.ylabel('Pressure Coefficient')
    plt.title('Pressure Coefficient Comparison')
    plt.xlim(-180+angle_of_attack, 180+angle_of_attack)
    plt.ylim(-3.5, 1.5)
    plt.legend()
    plt.show()

def plot_fields(panel_method, circle):
    X, Y = np.meshgrid(
        np.linspace(
            start = panel_method.surface.ro.x - 2 * circle.radius,
            stop = panel_method.surface.ro.x + 2 * circle.radius,
            num = 5
        ),
        np.linspace(
            start = panel_method.surface.ro.y - 2 * circle.radius,
            stop = panel_method.surface.ro.y + 2 * circle.radius,
            num = 50
        ),
        indexing='ij'
    )
    
    panel_method.display_contour(X, Y, BodyFixed_FrameOfReference=False)
    panel_method.display_velocity_field(X, Y, BodyFixed_FrameOfReference=False)
    panel_method.display_streamlines(X, Y, BodyFixed_FrameOfReference=False)
    
def test_BoundaryElementMethod():
    simulate_flow_around_a_2D_circular_object(
    radius=5,
    num_panels=10,
    center=(0, 0),
    velocity=1,
    angle_of_attack=0
)
