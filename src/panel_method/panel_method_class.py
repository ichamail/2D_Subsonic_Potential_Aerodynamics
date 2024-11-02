from __future__ import annotations
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
import numpy as np
from .mesh_class import AeroPanelMesh, SurfacePanelMesh
from ..myMath import Vector
from ..utilities import make_legend_arrow

class BoundaryElementMethod:
    
    def __init__(self, mesh:SurfacePanelMesh) -> None:
        self.mesh = mesh
        
        n = len(self.surface.face)
        
        self.B_ij = np.zeros((n, n))
        self.C_ij = np.zeros_like(self.B_ij)
        
        self.compute_surface_influence()
        
        self.V_inf = Vector(0, 0, 0)
                
        pass
    
    @property
    def surface(self):
        return self.mesh
    
    @property
    def panels(self):
        return self.surface.panel
    
    @property
    def A_ij(self):
        return self.C_ij
    
    @property
    def RHS(self):
        return - self.B_ij @ np.array(
            [panel.sigma for panel in self.surface.panel]
        )
    
    def set_V_inf(self, angle_of_attack, magnitude):
        
        self.set_BodyFixedFrame_orientation(0)
        self.set_BodyFixedFrame_angular_velocity(0)
        self.set_BodyFixedFrame_origin_velocity(0, 0)
        
        self.V_inf = Vector(
            magnitude * np.cos(np.deg2rad(angle_of_attack)),
            magnitude * np.sin(np.deg2rad(angle_of_attack)),
            0
        )
    
    def set_V_fs(self, angle_of_attack, magnitude):
        
        self.set_BodyFixedFrame_orientation(-angle_of_attack)
        self.set_BodyFixedFrame_angular_velocity(0)
        self.set_BodyFixedFrame_origin_velocity(Vo_x=-magnitude, Vo_y=0)
               
        self.V_inf = Vector(0, 0, 0)
    
    def V_fs_at_BodyFixedFrame_origin(self):
        return (self.V_inf - self.surface.Vo).transform(self.surface.A)    
    
    def V_fs_at(self, r_surface:Vector):
        
        V_fs_at_r_surface = (
            self.V_inf
            - (
                self.surface.Vo
                + self.surface.omega.cross(r_surface.transform(self.surface.A.T))
            )
        ).transform(self.surface.A)
        
        return V_fs_at_r_surface     
    
    def set_BodyFixedFrame_origin(self, xo, yo):
        self.surface.set_BodyFixedFrame_origin(xo, yo)
    
    def set_BodyFixedFrame_orientation(self, theta_z):
        self.surface.set_BodyFixedFrame_orientation(theta_z)
    
    def set_BodyFixedFrame_origin_velocity(self, Vo_x, Vo_y):
        self.surface.set_BodyFixedFrame_origin_velocity(Vo_x, Vo_y)
        
    def set_BodyFixedFrame_angular_velocity(self, omega_z):
        self.surface.set_BodyFixedFrame_angular_velocity(omega_z)      
            
    def compute_surface_influence(self):
        for panel_i in self.surface.panel:
            for panel_j in self.surface.panel:
                (
                    self.B_ij[panel_i.id][panel_j.id],
                    self.C_ij[panel_i.id][panel_j.id]
                ) = (
                    panel_j.unit_strength_induced_velocity_potential(
                        panel_i.r_cp
                    )
                )
    
    def compute_source_strengths(self):
        
        if self.surface.omega.norm() != 0:
            
            for panel in self.surface.panel:
                panel.sigma = - panel.e_n.dot(self.V_fs_at(panel.r_cp))
                
        else:
            V_fs = self.V_fs_at_BodyFixedFrame_origin()
            for panel in self.surface.panel:
                panel.sigma = - panel.e_n.dot(V_fs)
              
    def solve_linear_system(self):
             
        mu = np.linalg.solve(self.A_ij, self.RHS)
        
        for panel in self.surface.panel:
            panel.mu = mu[panel.id]
    
    def advance_solution(self):
       
       self.compute_source_strengths()
       
       self.solve_linear_system()
       
    def solve(self):
        self.advance_solution()
        self.compute_surface_velocity()
        self.compute_surface_pressure()
    
    def induced_velocity(self, r_p:Vector) -> Vector:
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: induced velocity vector of point P with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        v_induced = Vector(0, 0, 0)
        for panel in self.surface.panel:
            v_induced = v_induced + panel.induced_velocity(r_p)
        
        return v_induced       
    
    def induced_velocity_potential(self, r_p:Vector) -> float:
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            float: induced velocity potential at P
        """
        
        phi_induced = 0
        for panel in self.surface.panel:
            phi_induced = phi_induced + panel.induced_velocity_potential(r_p)
        
        return phi_induced 
    
    def velocity(self, r_p:Vector) -> Vector:
        
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: velocity vector of point P with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        if self.surface.omega.norm() != 0:
            V_fs = self.V_fs_at(r_p)                
        else:
            V_fs = self.V_fs_at_BodyFixedFrame_origin()               
              
        return V_fs + self.induced_velocity(r_p)
    
    def velocity_potential(self, r_p:Vector) -> float:
        
        """_summary_

        Args:
            r_p (Vector): position vector of point with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: velocity potential with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        if self.surface.omega.norm() != 0:
            print("error: unknow phi_infty when body is rotating")
            V_fs = self.V_fs_at(r_p)
            phi_fs = r_p.x * V_fs.x + r_p.y * V_fs.y                
        else:
            V_fs = self.V_fs_at_BodyFixedFrame_origin()
            phi_fs = r_p.x * V_fs.x + r_p.y * V_fs.y                
              
        return phi_fs + self.induced_velocity_potential(r_p)
       
    def compute_surface_velocity(self):
              
        for panel in self.surface.panel:
            
            # δεν λειτουργεί όπως και στην 3D έκδοση
            # panel.V = self.velocity(panel.r_cp)
            
            if self.surface.omega.norm() != 0:
                V_fs = self.V_fs_at(panel.r_cp)                
            else:
                V_fs = self.V_fs_at_BodyFixedFrame_origin()               
            
                            
            adjacent_panels = self.surface.give_adjacent_panels(panel)    
            
            if len(adjacent_panels) == 1:
                adjacent_panels = adjacent_panels + [
                    new_panel for new_panel in self.surface.give_adjacent_panels(
                        adjacent_panels[0]
                    ) if new_panel.id != panel.id 
                ]      
            
            panel.compute_surface_velocity(*adjacent_panels, V_fs)
          
    def compute_surface_pressure(self):
        
        if self.surface.omega.norm() != 0:
            
            for panel in self.surface.panel:
                V_fs = self.V_fs_at(panel.r_cp)
                panel.Cp = 1 - ( panel.V.norm()/V_fs.norm() )**2 
                             
        else:

            V_fs = self.V_fs_at_BodyFixedFrame_origin()
            for panel in self.surface.panel:
                panel.Cp = 1 - ( panel.V.norm()/V_fs.norm() )**2   
    
    def display_contour(self, X, Y, BodyFixed_FrameOfReference=True):
        nx, ny = X.shape
        phi = np.zeros_like(X, dtype=float)
        
        if BodyFixed_FrameOfReference:
            for i in range(nx):
                for j in range(ny):
                    if not self.surface.is_inside_polygon((X[i][j], Y[i][j])):
                        phi[i][j] = self.velocity_potential(
                            r_p=Vector(X[i][j], Y[i][j], 0)
                        )
        else:
            for i in range(nx):
                for j in range(ny):
                    r_p = (
                        Vector(X[i][j], Y[i][j], 0) - self.surface.ro
                    ).transform(self.surface.A)
                    if not self.surface.is_inside_polygon((r_p.x, r_p.y)):
                        phi[i][j] = self.velocity_potential(r_p)
                
        # plt.contour(
        #     X, Y, phi, levels=np.linspace(phi.min(), phi.max(), 20), zorder=1
        # )
        plt.contour(
            X, Y, phi, levels=20, zorder=1
        )
        
        if BodyFixed_FrameOfReference:
            plt.fill(
                self.surface.vertex[:,0], self.surface.vertex[:, 1], 'k', zorder=2
            )
        else:
            x = np.zeros(len(self.surface.vertex))
            y = np.zeros_like(x)
            
            for i, vertex in enumerate(self.surface.vertex):
                r = (
                    self.surface.ro 
                    + Vector(vertex[0], vertex[1], 0).transform(
                        self.surface.A.T
                    )
                )
                x[i] = r.x
                y[i] = r.y
                
            plt.fill(x, y, 'k', zorder=2)
        
        plt.axis('scaled')
        plt.show()
    
    def display_velocity_field(self, X, Y, BodyFixed_FrameOfReference=True):
        
        nx, ny = X.shape
        u = np.zeros_like(X, dtype=float)
        v = np.zeros_like(X, dtype=float)
        
        if BodyFixed_FrameOfReference:
            for i in range(nx):
                for j in range(ny):
                    if not self.surface.is_inside_polygon((X[i][j], Y[i][j])):
                        V = self.velocity(r_p=Vector(X[i][j], Y[i][j], 0))
                        u[i][j] += V.x
                        v[i][j] += V.y
        else:
            for i in range(nx):
                for j in range(ny):
                    r_p = (
                        Vector(X[i][j], Y[i][j], 0) - self.surface.ro
                    ).transform(self.surface.A)
                    if not self.surface.is_inside_polygon((r_p.x, r_p.y)):
                        V = self.velocity(r_p).transform(self.surface.A.T)
                        u[i][j] += V.x
                        v[i][j] += V.y
                        
        
        plt.quiver(X, Y, u, v, color='m', zorder=1)
        
        if BodyFixed_FrameOfReference:
            plt.fill(
                self.surface.vertex[:,0], self.surface.vertex[:, 1], 'k', zorder=2
            )
        else:
            x = np.zeros(len(self.surface.vertex))
            y = np.zeros_like(x)
            
            for i, vertex in enumerate(self.surface.vertex):
                r = (
                    self.surface.ro 
                    + Vector(vertex[0], vertex[1], 0).transform(
                        self.surface.A.T
                    )
                )
                x[i] = r.x
                y[i] = r.y
                
            plt.fill(x, y, 'k', zorder=2)
            
        plt.axis('scaled')
        plt.show()
                    
    def display_streamlines(self, X, Y, BodyFixed_FrameOfReference=True):
        nx, ny = X.shape
        u = np.zeros_like(X, dtype=float)
        v = np.zeros_like(X, dtype=float)
        
        extra_thickness_factor = 0.003 * abs(
            max(self.surface.vertex[:, 0]) - min(self.surface.vertex[:, 0])
        )
        self.surface.set_extra_thickness_layer(extra_thickness_factor)
        
        if BodyFixed_FrameOfReference:
            for i in range(nx):
                for j in range(ny):
                    if not self.surface.is_near_surface((X[i][j], Y[i][j])):
                        V = self.velocity(r_p=Vector(X[i][j], Y[i][j], 0))
                        u[i][j] += V.x
                        v[i][j] += V.y
        else:
            for i in range(nx):
                for j in range(ny):
                    r_p = (
                        Vector(X[i][j], Y[i][j], 0) - self.surface.ro
                    ).transform(self.surface.A)
                    if not self.surface.is_near_surface((r_p.x, r_p.y)):
                        V = self.velocity(r_p).transform(self.surface.A.T)
                        u[i][j] += V.x
                        v[i][j] += V.y
                
        start_points = np.column_stack([X[0, 3:-3:2], Y[0, 3:-3:2]])
        
        plt.streamplot(
            X.T, Y.T, u.T, v.T, start_points=start_points, linewidth=1.5,
            density = 50, arrowstyle='-', zorder=1
        )
        
                
        if BodyFixed_FrameOfReference:
            plt.fill(
                self.surface.vertex[:,0], self.surface.vertex[:, 1], 'k', zorder=2
            )
        else:
            x = np.zeros(len(self.surface.vertex))
            y = np.zeros_like(x)
            
            for i, vertex in enumerate(self.surface.vertex):
                r = (
                    self.surface.ro 
                    + Vector(vertex[0], vertex[1], 0).transform(
                        self.surface.A.T
                    )
                )
                x[i] = r.x
                y[i] = r.y
                
            plt.fill(x, y, 'k', zorder=2)
        
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Streamlines')   
        plt.axis('scaled')
        plt.show()
    

class PanelMethod(BoundaryElementMethod):
    
    def __init__(self, mesh:AeroPanelMesh) -> None:
        super().__init__(mesh)
        
        nw = len(self.wake.face)
        self.C_ij = np.pad(self.C_ij, ((0, 0), (0, nw)))
        
        self.steady_state = None
        
        
    @property
    def surface(self):
        return self.mesh.surface_mesh
    
    @property
    def wake(self):
        return self.mesh.wake_mesh
    
    @property
    def panels(self):
        return np.concatenate((self.surface.panel, self.wake.panel))
    
    @property
    def A_ij(self):
        
        Ns = len(self.surface.face)
        
        A_ij = np.array(self.C_ij[:, :Ns])
                
        if self.steady_state:
            
            for panel_i in self.surface.panel:   
                                    
                for panel_j in self.wake.panel:
            
                    A_ij[panel_i.id][self.surface.top_kutta_face_id] +=  self.C_ij[panel_i.id][Ns + panel_j.id]
                        
                    A_ij[panel_i.id][self.surface.bottom_kutta_face_id] -= self.C_ij[panel_i.id][Ns + panel_j.id]
            
        else:
            
            for panel_i in self.surface.panel:
                                
                A_ij[panel_i.id][self.surface.top_kutta_face_id] +=  self.C_ij[panel_i.id][Ns + self.wake.panel[-1].id]

                A_ij[panel_i.id][self.surface.bottom_kutta_face_id] -= self.C_ij[panel_i.id][Ns + self.wake.panel[-1].id]
        
        return A_ij
    
    @property
    def RHS(self):
        if self.steady_state:
            
            return super().RHS
        
        else:
            RHS = (
                - self.B_ij @ np.array(
                    [panel.sigma for panel in self.surface.panel]
                )
                - self.C_ij[:, len(self.surface.face):-1] @ np.array(
                    [panel.mu for panel in self.wake.panel[:-1]]
                )
            )
            return RHS 
    
    def set_V_inf(self, angle_of_attack, magnitude):
        
        self.set_BodyFixedFrame_orientation(0)
        self.set_BodyFixedFrame_angular_velocity(0)
        self.set_BodyFixedFrame_origin_velocity(0, 0)
        
        self.wake.set_BodyFixedFrame_orientation(angle_of_attack)
        self.V_inf = magnitude * Vector(1, 0, 0).transform(self.wake.A.T)
    
    def set_V_fs(self, angle_of_attack, magnitude):
        
        self.set_BodyFixedFrame_orientation(-angle_of_attack)
        self.set_BodyFixedFrame_angular_velocity(0)
        self.set_BodyFixedFrame_origin_velocity(Vo_x=-magnitude, Vo_y=0)
        
        if self.wake.surface_fixed:
            self.wake.set_BodyFixedFrame_orientation(angle_of_attack)
        else:
            self.wake.set_BodyFixedFrame_orientation(0)
            
        self.V_inf = Vector(0, 0, 0)
        
    def compute_wake_influence(self):
        
        Ns = len(self.surface.face)
        
        if self.wake.surface_fixed:
            
            for panel_i in self.surface.panel:
                
                for panel_j in self.wake.panel:
                    
                    r_cp = (panel_i.r_cp - self.wake.ro).transform(self.wake.A)
                    
                    self.C_ij[panel_i.id][Ns + panel_j.id] = (
                        panel_j.unit_strength_induced_velocity_potential(r_cp)
                    )
                                    
        else:
            
            for panel_i in self.surface.panel:
                
                for panel_j in self.wake.panel:
                    
                    r_cp = (
                        self.surface.ro 
                        + panel_i.r_cp.transform(self.surface.A.T)
                        - self.wake.ro
                    ).transform(self.wake.A)
                    
                    self.C_ij[panel_i.id][Ns + panel_j.id] = (
                        panel_j.unit_strength_induced_velocity_potential(r_cp)
                    )
                    
    def solve_linear_system(self):
               
        self.compute_wake_influence()
        
        super().solve_linear_system()
        
        mu_wake = (
                self.surface.panel[self.surface.top_kutta_face_id].mu
                - self.surface.panel[self.surface.bottom_kutta_face_id].mu
            )
        
        if self.steady_state:             
            for panel in self.wake.panel:
                panel.mu = mu_wake
        else:
           self.wake.panel[-1].mu = mu_wake
    
    def induced_velocity(self, r_p: Vector) -> Vector:
        v_induced = super().induced_velocity(r_p)
        
        if self.wake.surface_fixed:
            
            r_p = (r_p - self.wake.ro).transform(self.wake.A)
            for panel in self.wake.panel:
                v_induced = (
                    v_induced 
                    + panel.induced_velocity(r_p).transform(self.wake.A.T)
                )
        else:
            
            r_p = (
                self.surface.ro 
                + r_p.transform(self.surface.A.T)
                - self.wake.ro
            ).transform(self.wake.A)
            
            for panel in self.wake.panel:
                v_induced = (
                    v_induced 
                    + (panel.induced_velocity(r_p).transform(self.wake.A.T)
                    ).transform(self.surface.A)
                )            
                        
        return v_induced
    
    def induced_velocity_potential(self, r_p: Vector) -> float:
        phi_induced = super().induced_velocity_potential(r_p)
        
        if self.wake.surface_fixed:
            
            r_p = (r_p - self.wake.ro).transform(self.wake.A)
            for panel in self.wake.panel:
                phi_induced = phi_induced+ panel.induced_velocity_potential(r_p)
                
        else:
            
            r_p = (
                self.surface.ro 
                + r_p.transform(self.surface.A.T)
                - self.wake.ro
            ).transform(self.wake.A)
            
            for panel in self.wake.panel:
                phi_induced = phi_induced + (panel.induced_velocity_potential(r_p))
            
        return phi_induced   

    def solve(self, steady_state=True, iters=0):
        
        if steady_state:
            
            if iters==0:
                
                self.steady_state = steady_state
                return super().solve()
            
            else:
                
                return self.solve_iteratively(iters)
            
        else:
            
            V = self.V_fs_at(
                    Vector(
                        *self.surface.vertex[self.surface.kutta_vertex_id], 0
                    )
                ).norm()
            
            # length = (
            #     self.surface.panel[self.surface.top_kutta_face_id].length
            #     + self.surface.panel[self.surface.bottom_kutta_face_id].length
            # )/2
            
            length = max([panel.length for panel in self.surface.panel])
            
            dt = length/V 
            
                        
            self.solve_unsteady(dt, iters)
    
    def iter_wake(self):
        
        # positions vector's of wake's vertices with respect to body-fixed frame of reference f'
        if self.wake.surface_fixed:
                        
            r = np.array(
                [
                    self.wake.ro 
                    + Vector(*self.wake.vertex[i], 0).transform(self.wake.A.T)
                    for i in range(len(self.wake.vertex))
                ]
            )
            
        else:
            
            r = np.array(
                [
                    ( (self.wake.ro - self.surface.ro)
                      + Vector(*self.wake.vertex[i], 0).transform(self.wake.A.T)
                    ).transform(self.surface.A)
                    for i in range(len(self.wake.vertex))
                ]
            )
        
        
        ds = np.array(
            [(r[i+1] - r[i]).norm() for i in range(len(r)-1)]
        )
        
        i = 0
        v_s = (
            (self.surface.panel[self.surface.top_kutta_face_id].V
             + self.surface.panel[self.surface.bottom_kutta_face_id].V)/2
            + self.velocity(r[i+1])
        )/2
               
        e_s = v_s/v_s.norm()
        dr = ds[i] * e_s
        r[i+1] = r[i] + dr
        
        for i in range(1, len(r)-1):
            v_s = ( self.velocity(r[i]) + self.velocity(r[i+1]) )/2
            e_s = v_s/v_s.norm()
            dr = ds[i] * e_s
            r[i+1] = r[i] + dr
            
        
        # positions vector's of wake's vertices with respect to wake's frame of reference f''
        if self.wake.surface_fixed:
            
            for i in range(len(r)):
                r[i] = (r[i] - self.wake.ro).transform(self.wake.A)
        else:
            
            for i in range(len(r)):
                            
                r[i] = (
                    self.surface.ro + r[i].transform(self.surface.A.T)
                    - self.wake.ro
                ).transform(self.wake.A)
        
        for i in range(len(self.wake.vertex)):
            self.wake.vertex[i][0] = r[i].x
            self.wake.vertex[i][1] = r[i].y
        
        for face_id, face in enumerate(self.wake.face):
            for i, vertex_id in enumerate(face):
                self.wake.panel[face_id].r[i] = Vector(
                    *self.wake.vertex[vertex_id], 0
                )
       
    def solve_iteratively(self, iters=10):
        
        self.steady_state = True
        
        self.advance_solution()
        
        self.compute_surface_velocity()
        
        for i in range(iters):
            
            self.iter_wake()
            
            # self.mesh.display(
            #     BodyFixed_FrameOfReference=True, display_normals=True
            # )           
                        
            self.solve_linear_system()
            self.compute_surface_velocity()
        
        self.compute_surface_pressure()
            
    def set_WakeFixed_frame_origin_velocity(self):
        self.wake.set_BodyFixedFrame_origin_velocity(
            Vo_x=self.V_inf.x, Vo_y=self.V_inf.y
        )
    
    def solve_unsteady(self, dt:float, iters:int) -> None:
        
        self.steady_state = False
        
        self.set_WakeFixed_frame_origin_velocity()
                
        for i in range(iters):
            
            self.mesh.move_rigid_body(dt)
            
            # self.mesh.display(
            #     BodyFixed_FrameOfReference=False, display_normals=True
            # )
            
            self.C_ij = np.pad(self.C_ij, ((0, 0), (0, 1)))
                        
            self.advance_solution()
            
            self.roll_wake(dt)
            
            # self.mesh.display(
            #     BodyFixed_FrameOfReference=False, display_normals=True
            # )
            
            
        # self.mesh.display(
        #     BodyFixed_FrameOfReference=False, display_normals=True
        # )
        
        self.compute_surface_velocity()
        self.compute_surface_pressure()
    
    def roll_wake(self, dt):
        
        if self.wake.surface_fixed:
            print("wake must be fixed on a fluid particle in the undisturbed flow region")
                   
                  
        dr = np.array(
            [
                (
                    self.induced_velocity(
                        (
                            (self.wake.ro - self.surface.ro)
                            + Vector(*self.wake.vertex[i], 0).transform(
                                self.wake.A.T
                            )
                        ).transform(self.surface.A)
                    ).transform(self.surface.A.T)
                ).transform(self.wake.A) * dt
                for i in range(len(self.wake.vertex)-1)
            ]
        )       
                
        # for i in range(len(self.wake.vertex)-1):
        #     self.wake.move_vertex(vertex_id=i, dr=dr[i])     
                
        self.wake.move_vertices(dr)
        
    def AerodynamicForce(self, CharLength:float) -> Vector:
        CF = Vector(0, 0, 0)
        for panel in self.surface.panel:
            CF = CF + ( - panel.e_n * panel.Cp ) * (panel.length/CharLength)
        
        return CF
    
    def InducedDragForce(self, CharLength:float) -> Vector:
        
        CF = self.AerodynamicForce(CharLength)
        V_fs = self.V_fs_at_BodyFixedFrame_origin()
        e_fs = V_fs/V_fs.norm()
        
        CD = CF.dot(e_fs) * e_fs
        
        return CD
        
    def LiftForce(self, CharLength:float) -> Vector:
        
        CL = (
            self.AerodynamicForce(CharLength) 
            - self.InducedDragForce(CharLength)
        )
        
        return CL
    
    def Cm_about_point(self, r_p, CharLength:float) -> Vector:
        """_summary_
        
        Cm = Σ(r_i X CF_i) = Σ{(r_cp_i - r_p) X CF_i}
        
        Cm: moment coefficient about point p
        CF_i = n_i * (-Cp_i * A_i/A_ref)
        Cp_i: i-th's panel pressure coefficient
        A_i: i-th's panel area
        A_ref: Reference area
        r_cp_i: position vector of i-th's panel collocation point with respect to body-fixed frame of reference f'
        r_p: position vector of point P with respect to body-fixed frame of reference f'.
        r_cp_i = r_p + r_i
        

        Args:
            r_p (Vector): position vector of point with respect to body-fixed frame of reference f'. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference f' 
            
            CharLength (float): 2D analog for Reference Area of a 3D object

        Returns:
            Vector: moment coefficient vector Cm about point p
        """
        
        Cm = Vector(0, 0, 0)
        
        for panel in self.surface.panel:
            CF_i = (- panel.e_n * panel.Cp) * (panel.length/CharLength)
            r = panel.r_cp - r_p
            Cm = Cm + r.cross(CF_i)
        
        return Cm

    def Center_of_Pressure(self, CharLength:float) -> Vector:
        """_summary_
        
        a X b = c => b = (c X a)/(a*a) + ka, k: arbitary scalar constant
        
        (r_CoP - r_p) X F = Σ{r_i X Fi} = Σ{ (r_cp_i - r_p) X Fi} = Cm_about_p
        
        r_p = 0 e_x + 0 e_y + 0 e_z
        
        r_CoP X F = Σ{r_cp_i X Fi} = Cm_o  
        => F X r_CoP = - Σ{r_cp_i X Fi} = - Cm_o =>
        => r_CoP = (- Σ{r_cp_i X Fi} X F)/(F*F) + kF 
        => r_CoP = (F X Σ{r_cp_i X Fi})/(F*F) + kF
        => r_CoP = (CF X Cm_o)/(CF*CF) + k CF
        
        
        # r_cop * e_y = 0 =>  (CF X Cm_o)/(CF*CF) * e_y + kCF * e_y= 0 =>
        # k (F * e_y) = { (Cm_o X CF)/(CF*CF) } * e_y =>
        # k = { (Cm_o X CF)/(CF*CF) * e_y } / (CF * e_y)
        
        r_CoP: position vector of center of pressure, with respect to body-fixed frame of reference f'.
        F: Aerodynamic force = ΣFi
        r_cp_i: position vector of i-th's panel collocation point with respect to body-fixed frame of reference f' (position vector of the point where the force Fi act meassured from the body-fixed frame of reference f')
        r_p: position vector of point P with respect to body-fixed frame of reference f' (position vector of point P about which Cm is calculated).
        r_cp_i = r_p + r_i

        Args:
            CharLength (float): 2D analog for Reference Area of a 3D object

        Returns:
            Vector: position vector of center of pressure r_cop, with respect to body-fixed frame of reference f'. Coordinates of r_cop are expressed, also, with respect to the basis of body-fixed of reference f'
        """
        
        CF = self.AerodynamicForce(CharLength = CharLength)
        Cm_o = self.Cm_about_point(r_p=Vector(0, 0, 0), CharLength=CharLength)
        
        vec  = CF.cross(Cm_o)/(CF.norm()**2)
        k = - vec.y/CF.y
        
        r_cop = vec + k * CF
        
        return r_cop
        
    def display_forces(self) -> None:
        
        # ax, fig = self.mesh.plot(
        #     BodyFixed_FrameOfReference=False, display_normals=True
        # )
        
        ax, fig = self.surface.plot(
            BodyFixed_FrameOfReference=False, display_normals=False
        )

        CF = self.AerodynamicForce(CharLength = 1).transform(self.surface.A.T)
        r_cop = (
            self.surface.ro 
            + self.Center_of_Pressure(CharLength = 1).transform(self.surface.A.T)
        )
        
        ax.arrow(
            x=r_cop.x, y=r_cop.y, dx=CF.x, dy=CF.y,
            length_includes_head=True, head_width=0.025 * CF.norm(),
            color='m', label="$F$", zorder=3)
        
        ax.legend(
            handler_map={
                mpatches.FancyArrow: HandlerPatch(patch_func=make_legend_arrow)
            }
        )
        plt.show()
   