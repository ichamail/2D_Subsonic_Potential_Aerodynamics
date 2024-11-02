from src.panel_method import Circle, Airfoil


def plotCircle():
    
    Circle(
        name="circle", center=(0, 0), radius=1, num_points=11
    ).plot()
    
def plotAirfoil():
    
    name = "naca0012 sharp"
    chord = 2
    airfoil = Airfoil(name, chord)
    # airfoil.new_x_spacing(10)
    # airfoil.new_x_spacing2(5)
    # airfoil.new_x_spacing2(11, UpperLowerSpacing_equal=False)
    # print(airfoil.x_coords)
    # print(airfoil.y_coords)
    # airfoil.repanel(5+1, spacing="cosine")
    airfoil.repanel(5+1, spacing="denser at leading edge")
    airfoil.plot()
