import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from ..utilities import make_legend_arrow
from ..myMath import Vector


class Edge:
    def __init__(self, vertex_0, vertex_1, CCW=True) -> None:
        self._r_0 = Vector(*vertex_0, 0)
        self._r_1 = Vector(*vertex_1, 0)
        self.CCW = CCW
        self.compute_atrributes()
        
    @property
    def r_0(self):
        return self._r_0
    
    @r_0.setter
    def r_0(self, vector):
        self._r_0 = vector
        self.compute_atrributes()
    
    @property
    def r_1(self):
        return self._r_1
    
    @r_1.setter
    def r_1(self, vector):
        self._r_1 = vector
        self.compute_atrributes()
    
    def __getitem__(self, key):
        if key==0 or key==-1: return self.r_0
        if key==1 or key==-2: return self.r_1
        
    def __setitem__(self, key, vector:Vector) -> None:
        if key==0 or key==-1: self.r_0 = vector
        if key==1 or key==-2: self.r_1 = vector
    
    def compute_atrributes(self):
        self.e_t = self[1] - self[0]
        self.length = self.e_t.norm()
        self.e_t = self.e_t/self.length
        
        # e_z = Vector(0, 0, 1)
        # self.e_n = e_z.cross(self.e_t)
        # if self.CCW: self.e_n = -self.e_n
        
        e_z = Vector(0, 0, 1)
        if self.CCW:
            e_z = - e_z
            
        self.e_n = e_z.cross(self.e_t)
        
        
        self.A = np.array([[self.e_t.x, self.e_t.y, self.e_t.z],
                           [self.e_n.x, self.e_n.y, self.e_n.z],
                           [e_z.x, e_z.y, e_z.z]])
                
        self.r_cp = (self[0] + self[1])/2

        
class Panel:
    
    def __init__(self, vertex_0, vertex_1, CCW=True, id=-1) -> None:
        self.id = id
             
        # position vectors of panel's vertices
        self.edge = Edge(vertex_0, vertex_1, CCW)
                
        self.V = Vector(0, 0, 0) # surface velocity
        self.Cp = 0 # pressure coefficient
    
    @property
    def r(self):
        return self.edge
    
    @property
    def e_t(self):
        return  self.edge.e_t
    
    @property
    def e_n(self):
        return self.edge.e_n
    
    @property
    def A(self):
        return self.edge.A
    
    @property
    def r_cp(self):
        return self.edge.r_cp
    
    @property
    def length(self):
        return self.edge.length 
    
    def plot(self, PanelFixedFrame=False):
        
        fig = plt.figure()           
        ax = fig.add_subplot()
        
        if PanelFixedFrame:
            r_0 = (self.r[0] - self.r_cp).transform(self.A)
            r_1 = (self.r[1] - self.r_cp).transform(self.A)
            
            ax.plot([r_0.x, r_1.x], [r_0.y, r_1.y], "k", linewidth=2,zorder=1)
            
            # Panel-fixed frame of reference f'_j
            e_t = self.e_t.transform(self.A)
            e_n = self.e_n.transform(self.A)
            
            ax.arrow(
                0, 0, e_t.x, e_t.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{t_j}$", zorder=2
            )
            ax.arrow(
                0, 0, e_n.x, e_n.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{n_j}$", zorder=2
            )
            
            
            # Body-fixed frame of reference f'
            ro = -self.r_cp
            e_x = Vector(1, 0, 0).transform(self.A.T)
            e_y = Vector(0, 1, 0).transform(self.A.T)
            
            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label='$e_{x}$', zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label='$e_{y}$', zorder=2
            )
            
        else:
            
            ax.plot([self.r[0].x, self.r[1].x], [self.r[0].y, self.r[1].y], "k",linewidth=2, zorder=1)
            
            
            # Panel-fixed frame of reference f'_j       
            ax.arrow(
                self.r_cp.x, self.r_cp.y, self.e_t.x, self.e_t.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{t_j}$", zorder=2
            )
            ax.arrow(
                self.r_cp.x, self.r_cp.y, self.e_n.x, self.e_n.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{n_j}$", zorder=2
            )
            
            
            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0)
            e_y = Vector(0, 1, 0)
            
            ax.arrow(
                0, 0, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label='$e_{x}$', zorder=2
            )
            ax.arrow(
                0, 0, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label='$e_{y}$', zorder=2
            )
        
        ax.axis("equal")

        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        ax.legend(
                handler_map={
                    mpatches.FancyArrow: HandlerPatch(patch_func=make_legend_arrow)
                }
            )
        
        return ax, fig
    
    def display(self, PanelFixedFrame=False):
        
        ax, fig = self.plot(PanelFixedFrame)
        
        plt.show()

    def display_point_P(self, r_p:Vector, PanelFixedFrame=False):
        
        if PanelFixedFrame:
            r_p = (r_p - self.r_cp).transform(self.A) 
        
        ax, fig = self.plot(PanelFixedFrame)
        ax.scatter(r_p.x, r_p.y, c="r", edgecolors="k")
        plt.show()

           
class Source(Panel):
    
    def __init__(self, vertex_0, vertex_1, CCW=True, id=-1) -> None:
        super().__init__(vertex_0, vertex_1, CCW, id)
        
        self.sigma = 0  # source strength
    
    def unit_strength_induced_velocity_potential(self, r_p:Vector) -> float:
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        if r_p.norm() < 10**(-10):
           
            phi = self.length * np.log((self.length/2)**2)
                       
        elif r_p.y == 0 and r_0.x <= r_p.x <= r_1.x :
            
            if r_p.x == r_0.x or r_p.x == r_1.x:
                # r_p --> r_0 or r_1 and r_p inside panel => phi = -inf
                # r_p --> r_0 or r_1 and r_p outside panel => phi = inf 
                # phi = - np.inf
                phi = 0
            else:
                dr_0 = r_p - r_0
                dr_1 = r_p - r_1
                
                phi = (
                    dr_0.x * np.log(dr_0.norm()**2) 
                    - dr_1.x * np.log(dr_1.norm()**2)
                )
            
        else:
            
            dr_0 = r_p - r_0
            dr_1 = r_p - r_1
            
            phi = (
                dr_0.x * np.log(dr_0.norm()**2) 
                - dr_1.x * np.log(dr_1.norm()**2)
                + 2 * r_p.y * (
                    np.arctan2(r_p.y, dr_1.x) - np.arctan2(r_p.y, dr_0.x)
                )
            )
            
            phi = (
                (r_p.x - r_0.x) * np.log((r_0 - r_p).norm()**2)
                - (r_p.x - r_1.x) * np.log((r_1 - r_p).norm()**2)
                + 2 * r_p.y * ( np.arctan2(r_p.y, r_p.x - r_1.x) 
                            - np.arctan2(r_p.y, r_p.x - r_0.x) )
            )
                        
        return phi/(4 * np.pi)
                
    def induced_velocity_potential(self, r_p:Vector) -> float:
        return self.sigma * self.unit_strength_induced_velocity_potential(r_p)
    
    def induced_velocity(self, r_p:Vector) -> Vector:
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        
        
        if r_p.norm() == 0:
            
            # r_p.y --> 0^+  ==>  Vp.y --> np.pi
            # r_p.y --> 0^-  ==>  Vp.y --> - np.pi
            
            Vp = Vector(0, np.pi, 0)
        
        elif (r_p - r_0).norm() < self.length * 10**(-10) or (r_p - r_1).norm() < self.length * 10**(-10):
            # r_p --> r_0 or r_1   ==>  V_p.x --> inf
            # Vp = Vector(np.inf, np.pi, 0)
            # Vp = Vector(0, np.pi, 0)
            Vp = Vector(0, 0, 0)
            
        elif r_p.y == 0 and r_0.x <= r_p.x <= r_1.x :
            
            # r_p.y --> 0^+  ==>  V_p.y --> np.pi
            # r_p.y --> 0^-  ==>  V_p.y --> - np.pi
            
            if r_p.x == r_0.x or r_p.x == r_1.x:
                # r_p --> r_0 or r_1   ==>  V_p.x --> inf
                # Vp = Vector(np.inf, np.pi, 0)
                # Vp = Vector(0, np.pi, 0)
                Vp = Vector(0, 0, 0)
            
            else:
                
                Vp = Vector(
                    np.log((r_0 - r_p).norm()/(r_1 - r_p).norm()),
                    np.pi,
                    0
                )
                              
        else:
            
            Vp = Vector(
                np.log((r_0 - r_p).norm()/(r_1 - r_p).norm()),
                np.arctan2(r_p.y, r_p.x - r_1.x) 
                - np.arctan2(r_p.y, r_p.x - r_0.x),
                0
            )
        
        return ( self.sigma/(2*np.pi) * Vp ).transform(self.A.T)
                
    def induced_velocity(self, r_p:Vector) -> Vector:
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        norm_0 = (r_0 - r_p).norm()
        norm_1 = (r_1 - r_p).norm()
        
        if (
            norm_0 <= self.length * 10**(-12) 
            or norm_1 <= self.length * 10**(-12)
        ):
            
            Vx = 0
            
        else:
            
            Vx = np.log(norm_0/norm_1)
            
            
        if abs(r_p.y) <= self.length * 10**(-12):
            
            Vy = np.pi
            
        else:
            
            Vy = (
                np.arctan2(r_p.y, r_p.x - r_1.x) 
                - np.arctan2(r_p.y, r_p.x - r_0.x)
            )
            
        
        Vp = Vector(Vx, Vy, 0)
        
        return ( self.sigma/(2*np.pi) * Vp ).transform(self.A.T)

    
class Doublet(Panel):
    
    def __init__(self, vertex_0, vertex_1, CCW=True, id=-1) -> None:
        super().__init__(vertex_0, vertex_1, CCW, id)
        
        self.mu = 0  # doublet strength 
    
    def unit_strength_induced_velocity_potential(self, r_p:Vector) -> float:
                
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
              
        if r_p.norm() < 10**(-10):
            phi = - np.pi
            
        if r_p.y == 0 and r_0.x <= r_p.x <= r_1.x :
            # r_p.y --> 0^+  ==>  phi --> np.pi
            # r_p.y --> 0^-  ==>  phi --> - np.pi
            phi = - np.pi
            
        else:
            
            phi = (
                np.arctan2(r_p.y, r_p.x - r_1.x) 
                - np.arctan2(r_p.y, r_p.x - r_0.x)
            )
        
        return  phi/(2 * np.pi)
    
    def induced_velocity_potential(self, r_p:Vector) -> float:
        return self.mu * self.unit_strength_induced_velocity_potential(r_p)
    
    def induced_velocity(self, r_p:Vector) -> Vector:
        
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        if r_p.norm() == 0:
            
            Vp = Vector(0, -2, 0)
        
        elif (r_p - r_0).norm() < self.length * 10**(-10) or (r_p - r_1).norm() < self.length * 10**(-10):
            # r_p --> r_0 or r_1   ==>  V_p.x, V_p.y --> inf
            # Vp = Vector(0, np.inf, 0)
            Vp = Vector(0, 0, 0)
                                           
        elif r_p.y == 0 and r_0.x <= r_p.x <= r_1.x :
               
            if r_p.x == r_0.x or r_p.x == r_1.x:
                # r_p --> r_0 or r_1   ==>  V_p.x, V_p.y --> inf
                # Vp = Vector(0, np.inf, 0)
                Vp = Vector(0, 0, 0)
                
            else:
                
                dr_0 = r_p - r_0
                dr_1 = r_p - r_1
                Vp = Vector(0, 1/dr_1.x - 1/dr_0.x, 0)
                  
        else:
            dr_0 = r_p - r_0
            dr_1 = r_p - r_1
            norm_0 = dr_0.norm()**2
            norm_1 = dr_1.norm()**2
                
            Vp = Vector(
                - (dr_1.y/norm_1 - dr_0.y/norm_0),
                dr_1.x/norm_1 - dr_0.x/norm_0,
                0
            )
                        
        return ( self.mu/(2*np.pi) * Vp ).transform(self.A.T)

    def induced_velocity(self, r_p:Vector) -> Vector:
                
        # position vectors with respect to panel's frame of reference, expressed with unit vectors of panel's frame of refence
        r_p = (r_p - self.r_cp).transform(self.A) 
        r_0 = (self.r[0] - self.r_cp).transform(self.A)
        r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        if r_p.norm() <= 10**(-6):
            
            Vp = Vector(0, -2, 0)
            
        else:
            dr_0 = r_p - r_0
            dr_1 = r_p - r_1
            norm_0 = dr_0.norm()**2
            norm_1 = dr_1.norm()**2
            
            if norm_0 <= 10**(-6):
                
                Vp = Vector( - dr_1.y/norm_1 , dr_1.x/norm_1 , 0 )
                
            elif norm_1 <= 10**(-6):
                
                Vp = Vector( dr_0.y/norm_0 , - dr_0.x/norm_0 , 0 )    
                
            else:
                
                Vp = Vector(
                    - (dr_1.y/norm_1 - dr_0.y/norm_0),
                    dr_1.x/norm_1 - dr_0.x/norm_0,
                    0
                )
        
        return ( self.mu/(2*np.pi) * Vp ).transform(self.A.T)
        

class SurfacePanel(Source, Doublet):
    
    def unit_strength_induced_velocity_potential(self, r_p: Vector) -> float:
        return (Source.unit_strength_induced_velocity_potential(self, r_p),
                Doublet.unit_strength_induced_velocity_potential(self, r_p))
    
    def induced_velocity_potential(self, r_p: Vector) -> float:
        phi = self.unit_strength_induced_velocity_potential(r_p)
        return self.sigma * phi[0] + self.mu * phi[1]
            
    def induced_velocity(self, r_p: Vector) -> Vector:
        return (
            Source.induced_velocity(self, r_p) 
            + Doublet.induced_velocity(self, r_p)
        )
        
    def compute_surface_velocity(
        self, adjacent_panel_0:Panel, adjacent_panel_1:Panel, V_fs:Vector
    ):
        # r_0 = (self.r[0] - self.r_cp).transform(self.A)
        # r_1 = (self.r[1] - self.r_cp).transform(self.A)
        
        r_cp_0 = (adjacent_panel_0.r_cp - self.r_cp).transform(self.A)
        r_cp_1 = (adjacent_panel_1.r_cp - self.r_cp).transform(self.A)
        
        second_order_central_finite_difference_scheme = False
        second_order_forward_finite_difference_scheme = False
        second_order_backward_finite_difference_scheme = False
        
        if r_cp_0.x < 0 < r_cp_1.x:
            second_order_central_finite_difference_scheme = True
            # r_cp_prev = r_cp_0
            # r_cp_next = r_cp_1
            panel_prev = adjacent_panel_0
            panel_next = adjacent_panel_1
            
        elif r_cp_1.x < 0 < r_cp_0.x:
            second_order_central_finite_difference_scheme = True
            # r_cp_prev = r_cp_1
            # r_cp_next = r_cp_0
            panel_prev = adjacent_panel_1
            panel_next = adjacent_panel_0
        
        elif r_cp_0.x < r_cp_1.x < 0:
            second_order_backward_finite_difference_scheme = True
            panel_prev = adjacent_panel_1
            panel_prev_prev = adjacent_panel_0
            
        elif r_cp_1.x < r_cp_0.x < 0:
            second_order_backward_finite_difference_scheme = True
            panel_prev = adjacent_panel_0
            panel_prev_prev = adjacent_panel_1
            
        elif 0 < r_cp_0.x < r_cp_1.x:
            second_order_forward_finite_difference_scheme = True
            panel_next = adjacent_panel_0
            panel_next_next = adjacent_panel_1
        
        elif 0 < r_cp_1.x < r_cp_0.x:
            second_order_forward_finite_difference_scheme = True
            panel_next = adjacent_panel_1
            panel_next_next = adjacent_panel_0
        
        
        if second_order_central_finite_difference_scheme:
            # (du/dx)_i = ( u_i+1 - u_i-1 )/( 2*dx ) for constant length panels
            
            # s_next = r_1.norm() + (r_cp_next - r_1).norm()
            s_next = self.length/2 + panel_next.length/2
            
            # s_prev = r_0.norm() + (r_cp_prev - r_0).norm()
            s_prev = self.length/2 + panel_prev.length/2
            
            denom = s_prev + s_next
            num = (
                (panel_next.mu - self.mu) * s_prev/s_next
                - (panel_prev.mu - self.mu) * s_next/s_prev
            )
        
        elif second_order_forward_finite_difference_scheme:
            #(du/dx)_i = ( - 3*u_i + 4*u_i+1 - u_i+2 )/( 2*dx ) for constant length panels
            
            s_next = self.length/2 + panel_next.length/2
            s_next_next = (
                self.length/2 + panel_next.length + panel_next_next.length/2
            )
            
            denom = s_next_next - s_next
            num = (
                (self.mu - panel_next_next.mu) * s_next/s_next_next
                - (self.mu - panel_next.mu) * s_next_next/s_next  
            )
        
        elif second_order_backward_finite_difference_scheme:
            
            # (du/dx)_i = ( 3*u_i - 4*u_i-1 + u_i-2 )/( 2*dx ) for constant length panels
            
            s_prev = self.length/2 + panel_prev.length/2
            s_prev_prev = (
                self.length/2 + panel_prev.length + panel_prev_prev.length/2
            )
            
            denom = s_prev_prev - s_prev
            num = (
                (self.mu - panel_prev.mu) * s_prev_prev/s_prev
                - (self.mu - panel_prev_prev.mu) * s_prev/s_prev_prev
            )
        
        del_mu = num/denom
               
        self.V = (
            Vector(del_mu, self.sigma, 0).transform(self.A.T) + V_fs
        )

           
class WakePanel(Doublet):
    
    def unit_strength_induced_velocity_potential(self, r_p: Vector) -> float:
        return super().unit_strength_induced_velocity_potential(r_p)
    
    def induced_velocity_potential(self, r_p: Vector) -> float:
        return super().induced_velocity_potential(r_p)    
        
    def induced_velocity(self, r_p: Vector) -> Vector:
        return super().induced_velocity(r_p)
