from matplotlib import pyplot as plt
import numpy as np
from utilities import cosspace, DenserAtLeadingEdge, DenserAtTrailingEdge, is_inside_polygon
from scipy import interpolate
import sys


class Polygon:
    
    def __init__(self, name:str, coords, CCW_order=True) -> None:
        self.name = name
        self.coords = np.asarray(coords)
        
        if not CCW_order:
            self.invert_coords_order()
        
    def invert_coords_order(self):
        
        self.coords = np.array(
            [self.coords[-i] for i in range(len(self.coords))]
        )
      
    def is_inside_polygon(self, point:tuple) -> bool:
        
        return is_inside_polygon(
            [(self.coords[i, 0], self.coords[i, 1]) for i in range(len(self.coords))],
            point
        )
    
    def plot(self):
        plt.plot(self.coords[:, 0], self.coords[:, 1])
        plt.plot(
            self.coords[:, 0], self.coords[:, 1], 'ko', markerfacecolor='r'
        )
        plt.axis('scaled')
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.title(self.name)
        plt.show()
        

class Circle(Polygon):
    def __init__(
        self, name:str, center:tuple, radius:float, num_points:int=10
    ) -> None:
        
        self.radius = radius
        
        thetas = np.linspace(0, 2*np.pi, num_points)
        coords = np.array([
            [center[0] + radius * np.cos(theta),
             center[1] + radius * np.sin(theta)] for theta in thetas
        ])
        
        
        coords[-1] = coords[0]
            
        super().__init__(name, coords, CCW_order=True)


class Airfoil(Polygon):
    
    filePath="Airfoils/"
    
    def __init__(self, name:str, chord_length:float=1,
                 coords=None, num_points=100, CCW_order=True) -> None:
        
        if coords is None:
            self.name = name
            self.chord = chord_length
            self.get_from_data_base()
            self.coords = self.chord * self.coords
            if not CCW_order:
                self.invert_coords_order()
        else:
            super().__init__(name, coords, CCW_order)
            self.chord = chord_length
        
        self.repanel(num_points)
         
    
    @staticmethod
    def load_airfoil(filePath, fileName, header_lines=1):
        
        # Load the data from the text file
        fileName = filePath + fileName
        dataBuffer = np.loadtxt(fileName, delimiter=' ', skiprows=header_lines)
        
        # Extract data from the loaded dataBuffer array
        return dataBuffer
    
    def get_from_data_base(self):
        
        fileName = self.name + ".dat"
        try:
            self.coords = self.load_airfoil(self.filePath, fileName)
        except Exception as error:
            print("An error occurred:", type(error).__name__, "–", error)
            print("Airfoil doesn't exist in the database")
            sys.exit()
    
    def invert_coords_order(self):
        
        self.coords = np.array(
            [self.coords[-i] for i in range(len(self.coords))]
        )
        
    
    def close_trailing_edge(self):
        if self.coords[0] != self.coords[-1]:
            self.coords = np.vstack((self.coords, self.coords[0]))
    
    def give_suctionSide(self):
        index = np.where(self.coords[:, 0]==self.coords[:, 0].min())[0][0]
        return self.coords[0:index+1]
        
    
    def give_pressureSide(self):
        index = np.where(self.coords[:, 0]==self.coords[:, 0].min())[0][0]
        return self.coords[index:]
    
    
    def new_x_spacing(self, num_x_points):
        
        x, y = self.coords[:, 0], self.coords[:, 1]
          
        # Circle creation with diameter equal to airfoil chord
        x_max, x_min = max(x), min(x)
        R = (x_max - x_min)/2
        x_center = (x_max + x_min)/2
        theta = np.linspace(0, 2*np.pi, num_x_points+1)
        x_circle = x_center + R * np.cos(theta)
        
        # project circle points on x-axis
        x_project = np.copy(x_circle) # projections of x-cordiantes on airfoil
        y_project = np.empty_like(x_project)
        
        # compute y_project with interpolation
        j=0
        for i in range(num_x_points):
            while j < len(x)-1:
                if (x[j]<=x_project[i]<=x[j+1] or x[j+1]<=x_project[i]<=x[j]):
                    break
                else:
                    j = j+1
                
            # when break interpolate
            a = (y[j+1]-y[j])/(x[j+1]-x[j])
            b = y[j+1] - a * x[j+1]
            y_project[i] = a * x_project[i] + b
                
        y_project[num_x_points] = y_project[0]
        
        X, Y = x_project, y_project
        
        self.coords = np.column_stack([X, Y])
        
        
    def new_suctionSide_x_spacing(self, num_x_points):
        
        coords = self.give_suctionSide()
        x, y = coords[:, 0], coords[:, 1]
          
        # Circle creation with diameter equal to airfoil chord
        x_max, x_min = max(x), min(x)
        R = (x_max - x_min)/2
        x_center = (x_max + x_min)/2
        theta = np.linspace(0, np.pi, num_x_points+1)
        x_circle = x_center + R * np.cos(theta)
        
        # project circle points on x-axis
        x_project = np.copy(x_circle) # projections of x-cordiantes on airfoil
        y_project = np.empty_like(x_project)
        
        # compute y_project with interpolation
        j=0
        for i in range(num_x_points):
            while j < len(x)-1:
                if (x[j]<=x_project[i]<=x[j+1] or x[j+1]<=x_project[i]<=x[j]):
                    break
                else:
                    j = j+1
                
            # when break interpolate
            a = (y[j+1]-y[j])/(x[j+1]-x[j])
            b = y[j+1] - a * x[j+1]
            y_project[i] = a * x_project[i] + b
                
        y_project[num_x_points] = y_project[0]
        
        X, Y = x_project, y_project
               
        return np.column_stack([X, Y])
    
    def new_pressureSide_x_spacing(self, num_x_points):
        
        x, y = self.give_pressureSide()
          
        # Circle creation with diameter equal to airfoil chord
        x_max, x_min = max(x), min(x)
        R = (x_max - x_min)/2
        x_center = (x_max + x_min)/2
        theta = np.linspace(np.pi, 2*np.pi, num_x_points+1)
        x_circle = x_center + R * np.cos(theta)
        
        # project circle points on x-axis
        x_project = np.copy(x_circle) # projections of x-cordiantes on airfoil
        y_project = np.empty_like(x_project)
        
        # compute y_project with interpolation
        j=0
        for i in range(num_x_points):
            while j < len(x)-1:
                if (x[j]<=x_project[i]<=x[j+1] or x[j+1]<=x_project[i]<=x[j]):
                    break
                else:
                    j = j+1
                
            # when break interpolate
            a = (y[j+1]-y[j])/(x[j+1]-x[j])
            b = y[j+1] - a * x[j+1]
            y_project[i] = a * x_project[i] + b
                
        y_project[num_x_points] = y_project[0]
        
        X, Y = x_project, y_project
        
        return np.column_stack([X, Y])
        
    
    def new_x_spacing2(self, num_x_points, UpperLowerSpacing_equal=True):
        if UpperLowerSpacing_equal:
            SS_coords = self.new_suctionSide_x_spacing(num_x_points)
            PS_coords= self.new_pressureSide_x_spacing(num_x_points)
        else:
            quotient, remainder = np.divmod(num_x_points, 2)
            SS_coords= self.new_suctionSide_x_spacing(quotient)
            PS_coords = self.new_pressureSide_x_spacing(quotient+remainder)
        
        self.coords = np.vstack((SS_coords, PS_coords[1:])) 

    def repanel(self, n_points_per_side:int, spacing="cosine"):
        """
        improved version of new_x_spacing() and new_x_spacing2()
        """
        
        # upper and lower side coordinates
        SS_coords, PS_coords = self.give_suctionSide(), self.give_pressureSide()
        x_u, y_u = SS_coords[:, 0], SS_coords[:, 1]
        x_l, y_l = PS_coords[:, 0], PS_coords[:, 1]
        
        # distances between coordinates
        dr_u = np.sqrt( (x_u[:-1] - x_u[1:])**2 + (y_u[:-1] - y_u[1:])**2 )
        dr_l = np.sqrt( (x_l[:-1] - x_l[1:])**2 + (y_l[:-1] - y_l[1:])**2 )
        
        # distances from trailing edge
        dr_u = np.hstack((0, np.cumsum(dr_u)))
        dr_l = np.hstack((0, np.cumsum(dr_l)))
        
        # normalize
        dr_u = dr_u/dr_u[-1]
        dr_l = dr_l/dr_l[-1]
        
        dr = np.hstack((dr_u, 1 + dr_l[1:]))
        
        if spacing == "cosine":
            space = lambda  n_points_per_side: cosspace(0, 1, n_points_per_side) 
        elif spacing == "uniform":
            space = lambda n_points_per_side: np.linspace(0, 1, n_points_per_side)
        elif spacing == "denser at leading edge":
            space = lambda n_points_per_side: DenserAtLeadingEdge(0, 1, n_points_per_side, factor=1.4)
        elif spacing == "denser at trailing edge":
            space = lambda n_points_per_side: DenserAtTrailingEdge(0, 1, n_points_per_side, factor=1.4)
        else:
            space = lambda  n_points_per_side: cosspace(0, 1, n_points_per_side)
        
        # cosine-spaced list of points from 0 to 1
        x = space(n_points_per_side)
        # s = np.hstack((x, 1 + x[1:]))
        s = np.hstack((x, 2-x[-2::-1]))
        
        # Check that there are no duplicate points in the airfoil.
        if np.any(np.diff(dr))==0:
            raise ValueError(
                "This airfoil has a duplicated point (i.e. two adjacent points with the same (x, y) coordinates), so you can't repanel it!"
            )
        
        x_coords = interpolate.PchipInterpolator(dr, self.coords[:, 0])(s)
        y_coords = interpolate.PchipInterpolator(dr, self.coords[:, 1])(s)
        
        self.coords = np.column_stack([x_coords, y_coords])
    
        
    # def repanel(self, n_points_per_side:int, spacing="cosine"):
    #     """
    #     improved version of repanel. Needs some work
    #     """


    #     if spacing == "cosine":
    #         space = lambda start, stop, n_points_per_side: cosspace(start, stop, n_points_per_side) 
    #     elif spacing == "uniform":
    #         space = lambda start, stop, n_points_per_side: np.linspace(start, stop, n_points_per_side)
    #     # elif spacing == "denser at leading edge":
    #     #     space = lambda n_points_per_side: DenserAtLeadingEdge(n_points_per_side, factor=1.4)
    #     # elif spacing == "denser at trailing edge":
    #     #     space = lambda n_points_per_side: DenserAtTrailingEdge(n_points_per_side, factor=1.4)
    #     # else:
    #     #     space = lambda  n_points_per_side: cosspace(0, 1, n_points_per_side) 


    #     # upper and lower side coordinates
        
    #     old_upper_coordinates = self.give_suctionSide()
    #     old_lower_coordinates = self.give_pressureSide()

    #     # Find the streamwise distances between coordinates, assuming linear interpolation
    #     upper_distances_between_points = np.linalg.norm(np.diff(old_upper_coordinates, axis=0), axis=1)
    #     lower_distances_between_points = np.linalg.norm(np.diff(old_lower_coordinates, axis=0), axis=1)
    #     upper_distances_from_TE = np.concatenate(([0], np.cumsum(upper_distances_between_points)))
    #     lower_distances_from_LE = np.concatenate(([0], np.cumsum(lower_distances_between_points)))


    #     try:
    #         new_upper_coordinates = interpolate.CubicSpline(
    #             x=upper_distances_from_TE,
    #             y=old_upper_coordinates,
    #             axis=0,
    #             bc_type=(
    #                 (2, (0, 0)),
    #                 (1, (0, -1)),
    #             )
    #         )(space(0, upper_distances_from_TE[-1], n_points_per_side))

    #         new_lower_coordinates = interpolate.CubicSpline(
    #             x=lower_distances_from_LE,
    #             y=old_lower_coordinates,
    #             axis=0,
    #             bc_type=(
    #                 (1, (0, -1)),
    #                 (2, (0, 0)),
    #             )
    #         )(space(0, lower_distances_from_LE[-1], n_points_per_side))

    #     except ValueError as e:
    #         if not (
    #                 (np.all(np.diff(upper_distances_from_TE)) > 0) and
    #                 (np.all(np.diff(lower_distances_from_LE)) > 0)
    #         ):
    #             raise ValueError(
    #                 "It looks like your Airfoil has a duplicate point. Try removing the duplicate point and "
    #                 "re-running Airfoil.repanel()."
    #             )
    #         else:
    #             raise e

    #     self.coords = np.vstack(
    #         (new_upper_coordinates, new_lower_coordinates)
    #     )

        
    
if __name__=="__main__":
    
    name = "naca0012 sarp"
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
    
    
    Circle(
        name="circle", center=(0, 0), radius=1, num_points=11
    ).plot()
    
    