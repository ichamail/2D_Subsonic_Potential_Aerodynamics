from src.get_from_user import *
from tests.test_utilities import *
from tests.test_panel_method import *

def test():
    
    # # test spacing functions
    # test_spacing()
    
    # # test geometry class
    # plotCircle()
    # plotAirfoil()
    
    # # test panel class
    # test_panel()
    
    # # test mesh class
    # test_mesh()
    # test_aero_mesh()
    # test_aero_panel_mesh()
    # test_panel_mesh()
    # test_surface_mesh()
    # test_surface_panel_mesh()
    # test_wake_kinematics()
    # test_wake_mesh()
    # test_wake_panel_mesh()
    
    # # test panel method class
    # test_BoundaryElementMethod()
    # test_PanelMethod_SteadyState_rigidWake()
    # test_PanelMethod_SteadyState_iterativeWake()
    test_PanelMethod_Unsteady()
    
    pass
       
def main():
    
    # test()
    
    
    userInput = 0
    
    while userInput != 1 and userInput !=2:
        
        print(
            "\n1.Potential Flow around a Circle \n"
            + "2.Potential Flow around an Airfoil"
        )
                
        try:
            
            userInput = int(
                input(
                    "Select one from the above test cases by typing their corresponding number, and then press Enter:"
                )
            )
            
        except:
            
            print("\nplease insert a correct value")
    
    
    if userInput == 1:                 
        
        simulate_flow_around_a_2D_circular_object(
            radius = getLengthFromUser(name="Circle's radius", symbol="r"),
            center = getBodyFixedFrameFromUser(name="Circle's center"),
            velocity = getVelocityFromUser(),
            angle_of_attack = getAngleOfAttackFromUser(),
            num_panels = getNumOfPanelsFromUser(
                name="Suraface", symbol="Ns", min_panels=5, max_panels=200
            )
        )
    
    else:
        
        airfoil_name = getAirfoilNameFromUser(
            airfoil_list = ["naca0012 sharp"]
        )
        chord_length = getLengthFromUser(
            name="Airfoil's chord length", symbol="c"
        )
        leading_edge_location = getBodyFixedFrameFromUser(
            name="Airfoil's leading edge"
        )
        velocity = getVelocityFromUser()
        angle_of_attack = getAngleOfAttackFromUser()
        num_airfoil_panels = getNumOfPanelsFromUser(
            name="Surface", symbol="Ns", min_panels=5, max_panels=100
        )
        
        if getIfSteadyStateFromUser():
            
            wake_length_in_chords=getWakeLengthInChordsFromUser()
            
                        
            wake_relaxation_iters = getNumOfIterationsFromUser(
                type_of_iters="wake relaxation"
            )
            
            num_wake_panels = getNumOfPanelsFromUser(
                name="Wake", symbol="Nw",
                min_panels=wake_relaxation_iters, max_panels=100
            )
            
            print("\n\nSimulation Results:")
            
            simulate_steady_flow_around_an_Airfoil(
                airfoil_name = airfoil_name,
                chord_length = chord_length,
                leading_edge_location = leading_edge_location,
                velocity = velocity,
                angle_of_attack = angle_of_attack,
                num_airfoil_panels = num_airfoil_panels,
                wake_length_in_chords = wake_length_in_chords,
                num_wake_panels = num_wake_panels,
                kutta_vertex_id = 0,
                wake_relaxation_iters = wake_relaxation_iters 
            )
            
        
        else:
            
            print("\n\nSimulation Results:")
            
            simulate_unsteady_flow_around_an_Airfoil(
                airfoil_name=airfoil_name,
                chord_length=chord_length,
                leading_edge_location=leading_edge_location,
                velocity=velocity,
                angle_of_attack=angle_of_attack,
                num_airfoil_panels=num_airfoil_panels,
                kutta_vertex_id=0,
                num_time_steps = getNumOfIterationsFromUser(
                    type_of_iters="time"
                )
            )
        
        
    return 0

if __name__== "__main__":
    
    main()
