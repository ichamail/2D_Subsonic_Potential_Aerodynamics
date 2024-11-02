import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
import copy
from .panel_class import SurfacePanel, WakePanel, Panel
from ..utilities import make_legend_arrow, is_inside_polygon
from ..myMath import Vector


"""
F:X,Y,Z (e_X, e_Y, e_Z): inertial frame of reference
f:x,y,z (e_x, e_y, e_z): translating frame of reference
(e_X = e_x , e_Y = e_y, e_Z = e_z)

f':x',y',z' (e_x', e_y', e_z'): body-fixed frame of reference
    [e_x * e_x', e_y * e_x', e_z * e_x'
A =  e_x * e_y', e_y * e_y', e_z * e_y'
     e_x * e_z', e_y * e_z', e_z * e_z']

f":x",y",z" (e_x", e_y", e_z"): wake frame of reference

     [ e_x * e_x", e_y * e_x", e_z * e_x"
A_w =  e_x * e_y", e_y * e_y", e_z * e_y"    if orientation of f" is defined
       e_x * e_z", e_y * e_z", e_z * e_z" ]  relative to f'
     
     [ e_x * e_x", e_y * e_x", e_z * e_x"
A_w =  e_x * e_y", e_y * e_y", e_z * e_y"    if orientation of f" is defined  
       e_x * e_z", e_y * e_z", e_z * e_z" ]  relative to F


f_j': t_j,n_j,k_j (e_{t_j}, e_{n_j}, e_{k_j}): panel_j-fixed frame of reference 

     [ e_x' * e_t_j, e_y' * e_t_j, e_z' * e_t_j
A_j =  e_x' * e_n_j, e_y' * e_n_j, e_z' * e_n_j
       e_x' * e_k_j, e_y' * e_k_j, e_z' * e_k_j ]
"""

class Mesh:
    
    def __init__(self, vertex, face):
        
        
        self.ro = Vector(0, 0, 0) # position vector of mesh origin
        self.A = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]]) # orientation matrix A, A^T = R
        
        # velocity vector of body-fixed frame origin
        self.Vo = Vector(0, 0, 0)
        
        # body-fixed frame's angular velocity vector
        self.omega = Vector(0, 0, 0)
        
        self.vertex = np.copy(vertex)
        self.face = np.copy(face)
        
        self.adjacency_matrix = np.zeros((len(self.face), len(self.face)))
        self.set_adjacency_matrix()
        
    @property
    def vertex(self):
        return self._vertex
    
    @vertex.setter
    def vertex(self, arg):
        self._vertex = arg
    
    @property
    def face(self):
        return self._face
        
    @face.setter
    def face(self, arg):
        
        self._face = arg
        
    def set_BodyFixedFrame_origin(self, xo, yo):
        self.ro = Vector(xo, yo, 0)
        
    def set_BodyFixedFrame_orientation(self, theta_z):
        """
            A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
            A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) =>
            A(t+Δt)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx) A(t)]^T =>
            R(t+Δt) = A(t)^T Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t+Δt) = R(t) Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
            
            if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
            initial orientation of the body is such that f' coincides with f
            
            if A(t=0) =/= I then
            A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx) =>
            A(t=0)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx)]^T =>
            R(t=0) = Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t=0) = Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to finite rotations
            
            https://en.wikipedia.org/wiki/Infinitesimal_rotation_matrix
            https://en.wikipedia.org/wiki/Davenport_chained_rotations
        """
        theta_z = np.deg2rad(theta_z)
        Az = np.array([[np.cos(theta_z), np.sin(theta_z), 0],
                       [-np.sin(theta_z), np.cos(theta_z), 0],
                       [0, 0, 1]])
        self.A = np.copy(Az)
    
    def set_BodyFixedFrame_origin_velocity(self, Vo_x, Vo_y):
        self.Vo = Vector(Vo_x, Vo_y, 0)
        
    def set_BodyFixedFrame_angular_velocity(self, omega_z):
        self.omega = Vector(0, 0, np.deg2rad(omega_z))
        
    def move_BodyFixedFrame(self, dt):
        """
            A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
            A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) =>
            A(t+Δt)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx) A(t)]^T =>
            R(t+Δt) = A(t)^T Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t+Δt) = R(t) Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
            
            if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
            initial orientation of the body is such that f' coincides with f
            
            if A(t=0) =/= I then
            A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx) =>
            A(t=0)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx)]^T =>
            R(t=0) = Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t=0) = Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to finite rotations
            
            https://en.wikipedia.org/wiki/Infinitesimal_rotation_matrix
            https://en.wikipedia.org/wiki/Davenport_chained_rotations
        """
        
        self.ro = self.ro + self.Vo*dt
        dtheta =  self.omega*dt
        
        Az = np.array([[np.cos(dtheta.z), np.sin(dtheta.z), 0],
                       [-np.sin(dtheta.z), np.cos(dtheta.z), 0],
                       [0, 0, 1]])
        self.A = Az @ self.A
           
    def plot(self, BodyFixed_FrameOfReference=False):
        
        fig = plt.figure()           
        ax = fig.add_subplot()
        
        if BodyFixed_FrameOfReference:
            for face_id in range(len(self.face)):
                vertex_id = self.face[face_id][0]
                vertex_next_id = self.face[face_id][1]
                x,x_next = self.vertex[vertex_id][0], self.vertex[vertex_next_id][0]
                y,y_next = self.vertex[vertex_id][1], self.vertex[vertex_next_id][1]

                ax.plot([x, x_next], [y, y_next], "k", zorder=1)

            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0)  # e_x'
            e_y = Vector(0, 1, 0)  # e_y'

            ax.arrow(
                0, 0, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{x'}$", zorder=2
            )
            ax.arrow(
                0, 0, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{y'}$", zorder=2
            )

            # Inertial frame of reference F
            ro = -self.ro  # ro: r_oo' -> r_o'o = -roo'
            ro = ro.transform(self.A)
            e_x = Vector(1, 0, 0).transform(self.A)  # e_X
            e_y = Vector(0, 1, 0).transform(self.A)  # e_Y

            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label="$e_{x}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label="$e_{y}$", zorder=2
            )
            
        else:
            
            for face_id in range(len(self.face)):
                vertex_id = self.face[face_id][0]
                vertex_next_id = self.face[face_id][1]
                x, x_next = self.vertex[vertex_id][0], self.vertex[vertex_next_id][0]
                y,y_next = self.vertex[vertex_id][1], self.vertex[vertex_next_id][1]
                
                r = self.ro + Vector(x, y, 0).transform(self.A.T)
                r_next = self.ro + Vector(x_next, y_next, 0).transform(self.A.T)
        
                ax.plot([r.x, r_next.x], [r.y, r_next.y], "k", zorder=1)
            
            # Inertial frame of reference
            e_x = Vector(1, 0, 0)  # e_X
            e_y = Vector(0, 1, 0)  # e_Y

            ax.arrow(
                0, 0, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label="$e_{x}$", zorder=2
            )
            ax.arrow(
                0, 0, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label="$e_{y}$", zorder=2
            )
               
            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0).transform(self.A.T)  # e_x'
            e_y = Vector(0, 1, 0).transform(self.A.T)  # e_y'

            ro = self.ro
            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{x'}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{y'}$", zorder=2
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

    def display(self, BodyFixed_FrameOfReference=False):
        ax, fig = self.plot(BodyFixed_FrameOfReference)
        plt.show()
    
    def __copy__(self):
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj
    
    def copy(self):
        return self.__copy__()

    @staticmethod
    def do_intersect(face_i, face_j) -> bool:
        # intersections = 0
        # for vertex_id in face_i:
        #     if vertex_id in face_j:
        #         intersections = intersections + 1
                
        # if intersections >= 1 :
        #     return True
        # else:
        #     return False
        
        return sum(vertex_id in face_j for vertex_id in face_i)>=1  
        
    def set_adjacency_matrix(self):
        
        num_faces = len(self.face)
        self.adjacency_matrix = np.zeros((num_faces, num_faces))
        for i in range(num_faces):
            for j in range(num_faces):
                if i != j and self.do_intersect(self.face[i], self.face[j]):
                    self.adjacency_matrix[i][j] = 1
                                                
        # # should be faster but when meassured it isn't               
        # for i in range(num_faces):
        #     for j in range(num_faces):
        #         if self.adjacency_matrix[i][j] == 0:
        #             if i != j and self.do_intersect(self.face[i], self.face[j]):
        #                 self.adjacency_matrix[i][j] = 1
        #                 self.adjacency_matrix[j][i] = 1

        
class SurfaceMesh(Mesh):
    
    def __init__(self, vertex, face, kutta_vertex_id=None):
                
        super().__init__(vertex, face)
        
        self.kutta_vertex_id = kutta_vertex_id
        
        if kutta_vertex_id == None:
            self.top_kutta_face_id = None
            self.bottom_kutta_face_id = None          
        else:
            self.locate_kutta_faces_id()
            self.eliminate_trailing_edge_adjacency()
            
    
    def locate_kutta_faces_id(self):
        kutta_face = [
            id for id, face in enumerate(self.face) 
            if self.kutta_vertex_id in face
        ]
        self.top_kutta_face_id = kutta_face[0]
        self.bottom_kutta_face_id = kutta_face[1]
            
    def eliminate_trailing_edge_adjacency(self):
        self.adjacency_matrix[self.top_kutta_face_id][self.bottom_kutta_face_id] = 0
        self.adjacency_matrix[self.bottom_kutta_face_id][self.top_kutta_face_id] = 0

    def is_inside_polygon(self, point:tuple) -> bool:
        return is_inside_polygon(
            [
                (self.vertex[i, 0], self.vertex[i, 1])
                for i in range(len(self.vertex))
            ],
            point
        )
    
    def give_adjacent_faces(self, face_id):
        
        return [
            id for id in range(len(self.adjacency_matrix[face_id]))
            if self.adjacency_matrix[face_id][id] == 1
        ]

    
class WakeMesh(Mesh):

    def __init__(
        self, surface_mesh:SurfaceMesh, length, num_faces=1, surface_fixed=True
    ):
        

        self.surface_fixed = surface_fixed

        if self.surface_fixed:
            vertex = np.column_stack(
                (np.linspace(0, length, num_faces+1), np.zeros(num_faces+1))
            )
            face = np.array([[i, i+1] for i in range(num_faces)])
            
        else:
            vertex = np.column_stack(
                (np.linspace(length, 0, num_faces+1), np.zeros(num_faces+1))
            )
            
            face = np.array([[i+1, i] for i in range(num_faces)])
        
        super().__init__(vertex, face)
        
        
        if self.surface_fixed:
            """
            wake's frame of reference f'' is always fixed on surface's trailing edge
            
            wake's frame of refence orientation matrix A is defined relative to surface's body fixed frame of reference
            
            A = [e_x' * e_x", e_y' * e_x", e_z' * e_x"
                 e_x' * e_y", e_y' * e_y", e_z' * e_y"
                 e_x' * e_z", e_y' * e_z", e_z' * e_z"]
            
            r_o = x_o' e_x' + y_o' e_y' + z_o' e_z'  
            
            """
            self.set_BodyFixedFrame_origin(
                *surface_mesh.vertex[surface_mesh.kutta_vertex_id]
            )
        else:
            
            """
            wake's frame of reference f'' is not necessarily fixed on surface's trailing edge
            
            wake's frame of refence orientation matrix A is defined relative to inertial frame of reference F
            
            A = [e_X * e_x", e_Y * e_x", e_Z * e_x"
                 e_X * e_y", e_Y * e_y", e_Z * e_y"
                 e_X * e_z", e_Y * e_z", e_Z * e_z"]
            
            r_o = X_o e_X + Y_o e_Y + Z_o e_Z
            """
            
            ro = Vector(*surface_mesh.vertex[surface_mesh.kutta_vertex_id], 0)
            
            ro = surface_mesh.ro + ro.transform(surface_mesh.A.T)
            self.set_BodyFixedFrame_origin(ro.x, ro.y)
               
    def move_vertex(self, vertex_id:int, dr:Vector):
        self.vertex[vertex_id] = self.vertex[vertex_id] + (dr.x, dr.y)
    
    def add_vertex(self, x, y):
        self.vertex = np.vstack((self.vertex, [x, y]))
    
    def add_face(self, vertex_id, next_vertex_id):
        if len(self.face) == 0:
            self._face = np.array([[vertex_id, next_vertex_id]])
        else:
            self._face = np.vstack((self.face, [vertex_id, next_vertex_id]))

        self.append_adjacency_matrix()
          
    def append_adjacency_matrix(self):
        self.adjacency_matrix = np.pad(
            self.adjacency_matrix, pad_width=((0, 1), (0, 1)),
            mode="constant", constant_values=0
        )
        
        num_faces = len(self.face)
        i =  num_faces - 1
        for j in range(num_faces):
            if i != j and self.do_intersect(self.face[i], self.face[j]):
                self.adjacency_matrix[i][j] = 1
                self.adjacency_matrix[j][i] = 1
        
    def shed_wake(self, surface_mesh:Mesh):
        
        # orientation of f" must be defined relative to F 
        # (self.surface_fixed == False)
                     
        r_kutta = Vector(*surface_mesh.vertex[surface_mesh.kutta_vertex_id], 0)
        r_kutta = surface_mesh.ro + r_kutta.transform(surface_mesh.A.T)
        
        r_kutta = (r_kutta - self.ro).transform(self.A)
        
        self.add_vertex(r_kutta.x, r_kutta.y)
        self.add_face(len(self.vertex)-1, len(self.vertex)-2)      


class PanelMesh(Mesh):
    
    def __init__(self, vertex, face, CCW=True):
        self.CCW = CCW   
        super().__init__(vertex, face)
          
    @Mesh.face.setter
    def face(self, arg):
        
        self._face = arg
        self.panel = np.array(
            [
                Panel(
                    self.vertex[self.face[i][0]], self.vertex[self.face[i][1]],
                    CCW=self.CCW, id=i
                )
                for i in range(len(self.face))
            ]
        )
                
    @property
    def panel(self):
        return self._panel

    @panel.setter
    def panel(self, arg):
        self._panel = arg
                 
    def plot(self, BodyFixed_FrameOfReference=False, display_normals=False):
                
        if display_normals:
            
            fig = plt.figure()           
            ax = fig.add_subplot()
            
            if BodyFixed_FrameOfReference:
                
                # Inertial frame of reference F
                ro = -self.ro  # ro: r_oo' -> r_o'o = -roo'
                ro = ro.transform(self.A)
                e_x = Vector(1, 0, 0).transform(self.A)  # e_X
                e_y = Vector(0, 1, 0).transform(self.A)  # e_Y

                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='r', label="$e_{x}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='y', label="$e_{y}$", zorder=2
                )
                
                
                # Body-fixed frame of reference f'
                e_x = Vector(1, 0, 0)  # e_x'
                e_y = Vector(0, 1, 0)  # e_y'

                ax.arrow(
                    0, 0, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='b', label="$e_{x'}$", zorder=2
                )
                ax.arrow(
                    0, 0, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='g', label="$e_{y'}$", zorder=2
                )

                
                
                for panel in self.panel:
                    ax.plot(
                        [panel.r[0].x, panel.r[1].x], [panel.r[0].y, panel.r[1].y], "k", zorder=1
                    )
                    
                
                    ar_e_n = ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_n.x, panel.e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_t.x, panel.e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                
            else:
                
                # Inertial frame of reference
                e_x = Vector(1, 0, 0)  # e_X
                e_y = Vector(0, 1, 0)  # e_Y

                ax.arrow(
                    0, 0, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='r', label="$e_{x}$", zorder=2
                )
                ax.arrow(
                    0, 0, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='y', label="$e_{y}$", zorder=2
                )
                
                # Body-fixed frame of reference f'
                e_x = Vector(1, 0, 0).transform(self.A.T)  # e_x'
                e_y = Vector(0, 1, 0).transform(self.A.T)  # e_y'

                ro = self.ro
                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='b', label="$e_{x'}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='g', label="$e_{y'}$", zorder=2
                )
                
                for panel in self.panel:
                    
                    r_0 = self.ro + panel.r[0].transform(self.A.T)
                    r_1 = self.ro + panel.r[1].transform(self.A.T)
                    
                    ax.plot(
                        [r_0.x, r_1.x], [r_0.y, r_1.y], "k", zorder=1
                    )
                    
                                    
                    r_cp = self.ro + panel.r_cp.transform(self.A.T)
                    e_n = panel.e_n.transform(self.A.T)
                    e_t = panel.e_t.transform(self.A.T)
                    
                    ar_e_n = ax.arrow(
                        r_cp.x, r_cp.y, e_n.x, e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        r_cp.x, r_cp.y, e_t.x, e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )                   
                    
            
            ax.axis("equal")

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            if BodyFixed_FrameOfReference:
                ax.set_title("mesh displayed in body-fixed frame of reference f'")
            else:
                ax.set_title("mesh displayed in inertial frame of reference F")
                
            
            ar_e_t.set_label("$e_{t_j}$")
            ar_e_n.set_label("$e_{n_j}$")
                
            ax.legend(
                    handler_map={
                        mpatches.FancyArrow: HandlerPatch(patch_func=make_legend_arrow)
                    }
                )
            
            return ax, fig
        
        else:
            
            return super().plot(BodyFixed_FrameOfReference)
    
    def display(self, BodyFixed_FrameOfReference=False, display_normals=False):
        
        ax, fig = self.plot(BodyFixed_FrameOfReference, display_normals)
        
        plt.show()

    def copy(self):
        new_panel_mesh = super().copy()
        new_panel_mesh.panel = copy.deepcopy(self.panel)
        return new_panel_mesh


class SurfacePanelMesh(SurfaceMesh, PanelMesh):
    
    def __init__(self, vertex, face, CCW=True, kutta_vertex_id=None):
        
        self.CCW = CCW
        super().__init__(vertex, face, kutta_vertex_id)
        
        self.extra_thickness_layer_vertices = np.array([])

    @Mesh.face.setter
    def face(self, arg):
        
        self._face = arg
        self.panel = np.array(
            [
                SurfacePanel(
                    self.vertex[self.face[i][0]], self.vertex[self.face[i][1]],
                    CCW=self.CCW, id=i
                )
                for i in range(len(self.face))
            ]
        )

    def give_adjacent_panels(self, panel):
        return [self.panel[id] for id in super().give_adjacent_faces(panel.id)]

    
    def set_extra_thickness_layer(self, extra_thickness_factor:float):
        
        new_vertex = np.zeros_like(self.vertex)
        
        for vertex_id, vertex in enumerate(self.vertex):
            r = Vector(vertex[0], vertex[1], 0)
            panel = np.array([
                self.panel[face_id] for face_id in range(len(self.face))
                if vertex_id in self.face[face_id]
            ])
            dr = panel[0].e_n + panel[1].e_n
            dr = dr/dr.norm()
            r = r + extra_thickness_factor * dr
            new_vertex[vertex_id] = r.x, r.y
            
        self.extra_thickness_layer_vertices = new_vertex
    
    def is_near_surface(self, point:tuple) -> bool:
                
        if np.any(self.extra_thickness_layer_vertices):
            
            vertex = self.extra_thickness_layer_vertices
                                
            return is_inside_polygon(
                [(vertex[i, 0], vertex[i, 1])for i in range(len(vertex))], point
            )
            
        else:
            
            return is_inside_polygon(point)
            
    
    
class WakePanelMesh(PanelMesh, WakeMesh):
    
    def __init__(
        self, surface_mesh: SurfaceMesh, length, num_faces=1, surface_fixed=True
    ):   
        
        super(PanelMesh, self).__init__(surface_mesh, length, num_faces, surface_fixed)
        
                      
    @Mesh.face.setter
    def face(self, arg):
        
        self._face = arg
        self.panel = np.array(
            [
                WakePanel(
                    self.vertex[self.face[i][0]], self.vertex[self.face[i][1]],
                    CCW=False, id=i
                )
                for i in range(len(self.face))
            ]
        )     
    
    def move_vertex(self, vertex_id: int, dr: Vector):
        super().move_vertex(vertex_id, dr)
                
        for face_id_i in range(len(self.face)):
            for j, vertex_id_j in enumerate(self.face[face_id_i]):
                if vertex_id_j == vertex_id:
                    # self.panel[face_id_i].r[j] = self.panel[face_id_i].r[j] + dr
                    self.panel[face_id_i].r[j] = Vector(*self.vertex[vertex_id], 0)
    
    def move_vertices(self, dr:np.ndarray[any, Vector]):
        
        for i in range(len(dr)):
            super().move_vertex(vertex_id=i, dr=dr[i])
        
        for face_id, face in enumerate(self.face):
            for i, vertex_id in enumerate(face):
                self.panel[face_id].r[i] = Vector(*self.vertex[vertex_id], 0)
               
    def add_face(self, vertex_id, next_vertex_id, CCW:bool):
        super().add_face(vertex_id, next_vertex_id)
        self.panel = np.append(
            self.panel,
            WakePanel(self.vertex[self.face[-1][0]], self.vertex[self.face[-1][1]],
                    CCW=CCW, id=len(self.face)-1)
        )
    
    def shed_wake(self, surface_mesh: PanelMesh):
        # orientation of f" must be defined relative to F 
        # (self.surface_fixed == False)
                     
        r_kutta = Vector(*surface_mesh.vertex[surface_mesh.kutta_vertex_id], 0)
        r_kutta = surface_mesh.ro + r_kutta.transform(surface_mesh.A.T)
        
        r_kutta = (r_kutta - self.ro).transform(self.A)
        
        self.add_vertex(r_kutta.x, r_kutta.y)
        self.add_face(len(self.vertex)-1, len(self.vertex)-2, CCW=False)    
        
    
class AeroMesh:
    
    def __init__(self, surface_mesh, wake_mesh):
        self.surface_mesh = surface_mesh.copy()
        self.wake_mesh = wake_mesh.copy()
           
    def shed_wake(self):
        self.wake_mesh.shed_wake(self.surface_mesh)
    
    def move_rigid_body(self, dt):
        self.surface_mesh.move_BodyFixedFrame(dt)
        self.wake_mesh.move_BodyFixedFrame(dt) # wake_mesh.Vo = V_inf
        self.shed_wake()
        
    def plot(self, BodyFixed_FrameOfReference=True):
        fig = plt.figure()           
        ax = fig.add_subplot()
        
        if BodyFixed_FrameOfReference:
            
            # Inertial frame of reference F
            ro = -self.surface_mesh.ro  # ro: r_oo' -> r_o'o = -roo'
            ro = ro.transform(self.surface_mesh.A)
            e_x = Vector(1, 0, 0).transform(self.surface_mesh.A)  # e_X
            e_y = Vector(0, 1, 0).transform(self.surface_mesh.A)  # e_Y

            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label="$e_{x}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label="$e_{y}$", zorder=2
            )
            
            
            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0)  # e_x'
            e_y = Vector(0, 1, 0)  # e_y'

            ax.arrow(
                0, 0, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{x'}$", zorder=2
            )
            ax.arrow(
                0, 0, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{y'}$", zorder=2
            )
            
            
            # Wake's frame of reference f''
            if self.wake_mesh.surface_fixed:
                e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                ro = self.wake_mesh.ro
            else:
                e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                e_x = e_x.transform(self.surface_mesh.A)
                e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                e_y = e_y.transform(self.surface_mesh.A)
                
                ro = self.wake_mesh.ro - self.surface_mesh.ro
                ro = ro.transform(self.surface_mesh.A)
            
            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='m', label="$e_{x''}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='c', label="$e_{y''}$", zorder=2
            )

            
            
            # plot surface mesh
            for face_id in range(len(self.surface_mesh.face)):
                vertex_id = self.surface_mesh.face[face_id][0]
                vertex_next_id = self.surface_mesh.face[face_id][1]
                x = self.surface_mesh.vertex[vertex_id][0]
                y = self.surface_mesh.vertex[vertex_id][1]
                x_next = self.surface_mesh.vertex[vertex_next_id][0]
                y_next = self.surface_mesh.vertex[vertex_next_id][1]
                ax.plot([x, x_next], [y, y_next], "k", zorder=1)
            
            # plot wake mesh
            for face_id in range(len(self.wake_mesh.face)):
                vertex_id = self.wake_mesh.face[face_id][0]
                vertex_next_id = self.wake_mesh.face[face_id][1]
                x = self.wake_mesh.vertex[vertex_id][0]
                y = self.wake_mesh.vertex[vertex_id][1]
                x_next = self.wake_mesh.vertex[vertex_next_id][0]
                y_next = self.wake_mesh.vertex[vertex_next_id][1]
                
                if self.wake_mesh.surface_fixed:
                    r = self.wake_mesh.ro + Vector(x, y, 0).transform(self.wake_mesh.A.T)
                    r_next = self.wake_mesh.ro + Vector(x_next, y_next, 0).transform(self.wake_mesh.A.T)
                else:
                    
                    ro = self.wake_mesh.ro - self.surface_mesh.ro
                    
                    r = ro + Vector(x, y, 0).transform(
                        self.wake_mesh.A.T
                    ) 
                    
                    r_next = ro + Vector(x_next, y_next, 0).transform(
                        self.wake_mesh.A.T
                    ) 
                    
                    r = r.transform(self.surface_mesh.A)
                    r_next = r_next.transform(self.surface_mesh.A)
                                        
                ax.plot(
                    [r.x, r_next.x], [r.y, r_next.y], "cornflowerblue", zorder=1
                )
                
        else:
            
            # Inertial frame of reference F
            e_x = Vector(1, 0, 0)  # e_X
            e_y = Vector(0, 1, 0)  # e_Y

            ax.arrow(
                0, 0, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='r', label="$e_{x}$", zorder=2
            )
            ax.arrow(
                0, 0, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='y', label="$e_{y}$", zorder=2
            )
                
            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0).transform(self.surface_mesh.A.T)  # e_x'
            e_y = Vector(0, 1, 0).transform(self.surface_mesh.A.T)  # e_y'

            ro = self.surface_mesh.ro
            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='b', label="$e_{x'}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='g', label="$e_{y'}$", zorder=2
            )

            # Wake's frame of reference f''
            
            if self.wake_mesh.surface_fixed:
                ro = (
                    self.surface_mesh.ro
                    + self.wake_mesh.ro.transform(self.surface_mesh.A.T)
                )
                e_x = Vector(1, 0, 0).transform(
                        self.surface_mesh.A.T @ self.wake_mesh.A.T
                    )  # e_x"
                e_y = Vector(0, 1, 0).transform(
                        self.surface_mesh.A.T @ self.wake_mesh.A.T
                    )  # e_y"
                
            else:
                ro = self.wake_mesh.ro
                e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                
            ax.arrow(
                ro.x, ro.y, e_x.x, e_x.y,
                length_includes_head=True, head_width=0.05,
                color='m', label="$e_{x''}$", zorder=2
            )
            ax.arrow(
                ro.x, ro.y, e_y.x, e_y.y,
                length_includes_head=True, head_width=0.05,
                color='c', label="$e_{y''}$", zorder=2
            )
            
            
            # plot surface mesh
            for face_id in range(len(self.surface_mesh.face)):
                vertex_id = self.surface_mesh.face[face_id][0]
                vertex_next_id = self.surface_mesh.face[face_id][1]
                x = self.surface_mesh.vertex[vertex_id][0]
                x_next = self.surface_mesh.vertex[vertex_next_id][0]
                y = self.surface_mesh.vertex[vertex_id][1]
                y_next = self.surface_mesh.vertex[vertex_next_id][1]
                
                r = (
                    self.surface_mesh.ro 
                    + Vector(x, y, 0).transform(self.surface_mesh.A.T)
                )
                
                r_next = (
                    self.surface_mesh.ro 
                    + Vector(x_next, y_next, 0).transform(self.surface_mesh.A.T)
                )
        
                ax.plot([r.x, r_next.x], [r.y, r_next.y], "k", zorder=1)
            
            # plot wake mesh
            for face_id in range(len(self.wake_mesh.face)):
                vertex_id = self.wake_mesh.face[face_id][0]
                vertex_next_id = self.wake_mesh.face[face_id][1]
                x = self.wake_mesh.vertex[vertex_id][0]
                x_next = self.wake_mesh.vertex[vertex_next_id][0]
                y = self.wake_mesh.vertex[vertex_id][1]
                y_next = self.wake_mesh.vertex[vertex_next_id][1]
                
                
                if self.wake_mesh.surface_fixed:
                    
                    # position vector of wake vertex with respect to origin of body-fixed frame of reference f'               
                    r = (
                        self.wake_mesh.ro 
                        + Vector(x, y, 0).transform(self.wake_mesh.A.T)
                    )
                    
                    # position vector of wake vertex with respect to origin of inertial frame of reference F
                    r = (
                        self.surface_mesh.ro 
                        + r.transform(self.surface_mesh.A.T)
                    )
                
                    r_next = (
                        self.wake_mesh.ro 
                        + Vector(x_next, y_next, 0).transform(
                                                        self.wake_mesh.A.T
                                                    )
                    )
                    
                    r_next = (
                        self.surface_mesh.ro
                        + r_next.transform(self.surface_mesh.A.T)
                    )
                
                else:
                    # position vector of wake vertex with respect to origin of inertial frame of reference F
                    r = (
                        self.wake_mesh.ro 
                        + Vector(x, y, 0).transform(self.wake_mesh.A.T)
                    )
                    
                    r_next = (
                        self.wake_mesh.ro 
                        + Vector(x_next, y_next, 0).transform(
                                                        self.wake_mesh.A.T
                                                    )
                    )
        
                ax.plot(
                    [r.x, r_next.x], [r.y, r_next.y], "cornflowerblue", zorder=1
                )
            
            
        ax.axis("equal")

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        
        if BodyFixed_FrameOfReference:
            ax.set_title("mesh displayed in body-fixed frame of reference f'")
        else:
            ax.set_title("mesh displayed in inertial frame of reference F")
        
        ax.legend(
                handler_map={
                    mpatches.FancyArrow: HandlerPatch(patch_func=make_legend_arrow)
                }
            )
                
        return ax, fig
        
    def display(self, BodyFixed_FrameOfReference=True):
        ax, fig = self.plot(BodyFixed_FrameOfReference)
        plt.show()
    
class AeroPanelMesh(AeroMesh): 
    
    def plot(self, BodyFixed_FrameOfReference=True, display_normals=False):
        
        if display_normals:
            fig = plt.figure()           
            ax = fig.add_subplot()
            
            if BodyFixed_FrameOfReference:
                
                # Inertial frame of reference F
                ro = -self.surface_mesh.ro  # ro: r_oo' -> r_o'o = -roo'
                ro = ro.transform(self.surface_mesh.A)
                e_x = Vector(1, 0, 0).transform(self.surface_mesh.A)  # e_X
                e_y = Vector(0, 1, 0).transform(self.surface_mesh.A)  # e_Y

                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='r', label="$e_{x}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='y', label="$e_{y}$", zorder=2
                )
                
                # Body-fixed frame of reference f'
                e_x = Vector(1, 0, 0)  # e_x'
                e_y = Vector(0, 1, 0)  # e_y'

                ax.arrow(
                    0, 0, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='b', label="$e_{x'}$", zorder=2
                )
                ax.arrow(
                    0, 0, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='g', label="$e_{y'}$", zorder=2
                )
                
                
                # Wake's frame of reference f''
                if self.wake_mesh.surface_fixed:
                    e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                    e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                    ro = self.wake_mesh.ro
                else:
                    e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                    e_x = e_x.transform(self.surface_mesh.A)
                    e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                    e_y = e_y.transform(self.surface_mesh.A)
                    
                    ro = self.wake_mesh.ro - self.surface_mesh.ro
                    ro = ro.transform(self.surface_mesh.A)
                
                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='tab:blue', label="$e_{x''}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='tab:orange', label="$e_{y''}$", zorder=2
                )

                
                # plot surface mesh
                for panel in self.surface_mesh.panel:
                    ax.plot(
                        [panel.r[0].x, panel.r[1].x],
                        [panel.r[0].y, panel.r[1].y], "k", zorder=1
                    )
                
                    ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_n.x, panel.e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_t.x, panel.e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                    
                    ar_e_n = ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_n.x, panel.e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        panel.r_cp.x, panel.r_cp.y, panel.e_t.x, panel.e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                
                # plot wake mesh
                for panel in self.wake_mesh.panel:
                    
                    if self.wake_mesh.surface_fixed:
                        r = self.wake_mesh.ro + panel.r[0].transform(self.wake_mesh.A.T)
                        r_next = self.wake_mesh.ro + panel.r[1].transform(self.wake_mesh.A.T)
                        
                        r_cp = self.wake_mesh.ro + panel.r_cp.transform(self.wake_mesh.A.T)
                        
                        e_t = panel.e_t.transform(self.wake_mesh.A.T)
                        e_n = panel.e_n.transform(self.wake_mesh.A.T)
                        
                        
                    else:
                        
                        ro = self.wake_mesh.ro - self.surface_mesh.ro
                        
                        r = ro + panel.r[0].transform(
                            self.wake_mesh.A.T
                        ) 
                        
                        r_next = ro + panel.r[1].transform(
                            self.wake_mesh.A.T
                        ) 
                        
                        r_cp = ro + panel.r_cp.transform(self.wake_mesh.A.T)
                        e_t = panel.e_t.transform(self.wake_mesh.A.T)
                        e_n = panel.e_n.transform(self.wake_mesh.A.T)
                        
                        
                        r = r.transform(self.surface_mesh.A)
                        r_next = r_next.transform(self.surface_mesh.A)
                        r_cp = r_cp.transform(self.surface_mesh.A)
                        e_t = e_t.transform(self.surface_mesh.A)
                        e_n = e_n.transform(self.surface_mesh.A)
                                            
                    ax.plot([r.x, r_next.x], [r.y, r_next.y], "k", zorder=1)
                    
                    ar_e_n = ax.arrow(
                        r_cp.x, r_cp.y, e_n.x, e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        r_cp.x, r_cp.y, e_t.x, e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                    
                
            else:
                
                # Inertial frame of reference F
                e_x = Vector(1, 0, 0)  # e_X
                e_y = Vector(0, 1, 0)  # e_Y

                ax.arrow(
                    0, 0, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='r', label="$e_{x}$", zorder=2
                )
                ax.arrow(
                    0, 0, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='y', label="$e_{y}$", zorder=2
                )
                    
                # Body-fixed frame of reference f'
                e_x = Vector(1, 0, 0).transform(self.surface_mesh.A.T)  # e_x'
                e_y = Vector(0, 1, 0).transform(self.surface_mesh.A.T)  # e_y'

                ro = self.surface_mesh.ro
                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='b', label="$e_{x'}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='g', label="$e_{y'}$", zorder=2
                )

                # Wake's frame of reference f''
                
                if self.wake_mesh.surface_fixed:
                    ro = (
                        self.surface_mesh.ro
                        + self.wake_mesh.ro.transform(self.surface_mesh.A.T)
                    )
                    e_x = Vector(1, 0, 0).transform(
                            self.surface_mesh.A.T @ self.wake_mesh.A.T
                        )  # e_x"
                    e_y = Vector(0, 1, 0).transform(
                            self.surface_mesh.A.T @ self.wake_mesh.A.T
                        )  # e_y"
                    
                else:
                    ro = self.wake_mesh.ro
                    e_x = Vector(1, 0, 0).transform(self.wake_mesh.A.T)  # e_x"
                    e_y = Vector(0, 1, 0).transform(self.wake_mesh.A.T)  # e_y"
                    
                ax.arrow(
                    ro.x, ro.y, e_x.x, e_x.y,
                    length_includes_head=True, head_width=0.05,
                    color='tab:blue', label="$e_{x''}$", zorder=2
                )
                ax.arrow(
                    ro.x, ro.y, e_y.x, e_y.y,
                    length_includes_head=True, head_width=0.05,
                    color='tab:orange', label="$e_{y''}$", zorder=2
                )
                
                
                # plot surface mesh
                for panel in self.surface_mesh.panel:
                    
                    r = (
                        self.surface_mesh.ro 
                        + panel.r[0].transform(self.surface_mesh.A.T)
                    )
                    
                    r_next = (
                        self.surface_mesh.ro 
                        + panel.r[1].transform(self.surface_mesh.A.T)
                    )
            
                                        
                    r_cp = (
                        self.surface_mesh.ro 
                        + panel.r_cp.transform(self.surface_mesh.A.T)
                    )
                    
                    e_t =  panel.e_t.transform(self.surface_mesh.A.T)
                    
                    e_n = panel.e_n.transform(self.surface_mesh.A.T)
                    
                    ax.plot([r.x, r_next.x], [r.y, r_next.y], "k", zorder=1)
                    
                    ar_e_n = ax.arrow(
                        r_cp.x, r_cp.y, e_n.x, e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        r_cp.x, r_cp.y, e_t.x, e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                
                # plot wake mesh
                for panel in self.wake_mesh.panel:
                                       
                    if self.wake_mesh.surface_fixed:
                        
                        # position vector of wake vertex with respect to origin of body-fixed frame of reference f'               
                        r = (
                            self.wake_mesh.ro 
                            + panel.r[0].transform(self.wake_mesh.A.T)
                        )
                        
                        # position vector of wake vertex with respect to origin of inertial frame of reference F
                        r = (
                            self.surface_mesh.ro 
                            + r.transform(self.surface_mesh.A.T)
                        )
                    
                        r_next = (
                            self.wake_mesh.ro 
                            + panel.r[1].transform(self.wake_mesh.A.T)
                        )
                        
                        r_next = (
                            self.surface_mesh.ro
                            + r_next.transform(self.surface_mesh.A.T)
                        )
                        

                        r_cp = (
                            self.wake_mesh.ro
                            + panel.r_cp.transform(self.wake_mesh.A.T)
                        )
                        
                        r_cp = (
                            self.surface_mesh.ro
                            + r_cp.transform(self.surface_mesh.A.T)
                        )
                        
                        e_t = panel.e_t.transform(self.wake_mesh.A.T)
                        
                        e_t = e_t.transform(self.surface_mesh.A.T)
                        
                        
                        e_n = panel.e_n.transform(self.wake_mesh.A.T)
                        
                        e_n = e_n.transform(self.surface_mesh.A.T)
                        
                    else:
                        # position vector of wake vertex with respect to origin of inertial frame of reference F
                        r = (
                            self.wake_mesh.ro 
                            + panel.r[0].transform(self.wake_mesh.A.T)
                        )
                        
                        r_next = (
                            self.wake_mesh.ro 
                            + panel.r[1].transform(self.wake_mesh.A.T)
                        )
                        
                        r_cp = (
                            self.wake_mesh.ro
                            + panel.r_cp.transform(self.wake_mesh.A.T)
                        )
                        
                        e_t = panel.e_t.transform(self.wake_mesh.A.T)
                        
                        e_n = panel.e_n.transform(self.wake_mesh.A.T)
            
                    ax.plot([r.x, r_next.x], [r.y, r_next.y], "k", zorder=1)
                    
                    ar_e_n = ax.arrow(
                        r_cp.x, r_cp.y, e_n.x, e_n.y,
                        length_includes_head=True, head_width=0.05,
                        color='m', zorder=2
                    )
                    ar_e_t = ax.arrow(
                        r_cp.x, r_cp.y, e_t.x, e_t.y,
                        length_includes_head=True, head_width=0.05,
                        color='c', zorder=2
                    )
                    
                
            ax.axis("equal")

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            if BodyFixed_FrameOfReference:
                ax.set_title("mesh displayed in body-fixed frame of reference f'")
            else:
                ax.set_title("mesh displayed in inertial frame of reference F")
            
            
            ar_e_t.set_label("$e_{t_j}$")
            ar_e_n.set_label("$e_{n_j}$")
            
            ax.legend(
                    handler_map={
                        mpatches.FancyArrow: HandlerPatch(patch_func=make_legend_arrow)
                    }
                )
                    
            return ax, fig 
  
        else:
            return super().plot(BodyFixed_FrameOfReference)
    
    def display(self, BodyFixed_FrameOfReference=True, display_normals=False):
        ax, fig = self.plot(BodyFixed_FrameOfReference, display_normals)
        plt.show()    
