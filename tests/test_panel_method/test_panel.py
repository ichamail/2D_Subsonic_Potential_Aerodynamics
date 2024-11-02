import numpy as np
from matplotlib import pyplot as plt
from src.myMath import Vector
from src.panel_method import Panel, Doublet, Source


def test_panel(
    coords:tuple[tuple[float, float]]=((3, 1), (1, 2)),
    CCW:bool=True,
    r_p:Vector=Vector(4, 3, 0)
)-> None:
    
    panel = Panel(coords[0], coords[1], CCW)
    
    doublet_panel = Doublet(coords[0], coords[1], CCW)
    doublet_panel.mu = 1
    
    source_panel = Source(coords[0], coords[1], CCW)
    source_panel.sigma = 1
    
    panel.display_point_P(r_p, PanelFixedFrame=False)
    panel.display_point_P(r_p, PanelFixedFrame=True)
    
    
    x = np.linspace(panel.r_cp.x - 4 * panel.length,
                    panel.r_cp.x + 4 * panel.length, 100)
    y = np.linspace(panel.r_cp.y - 4 * panel.length,
                    panel.r_cp.y + 4 * panel.length, 100)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    u_src = np.zeros_like(X, dtype=float)
    v_src = np.zeros_like(X, dtype=float)
    phi_src = np.zeros_like(X, dtype=float)
    
    u_dblt = np.zeros_like(X, dtype=float)
    v_dblt = np.zeros_like(X, dtype=float)
    phi_dblt = np.zeros_like(X, dtype=float)
    
    
    for i in range(len(x)):
        for j in range(len(y)):
            
            phi_src[i][j] = source_panel.induced_velocity_potential(
                Vector(X[i][j], Y[i][j], 0)
            )
            
            phi_dblt[i][j] = doublet_panel.induced_velocity_potential(
                Vector(X[i][j], Y[i][j], 0)
            )
            
            V_src = source_panel.induced_velocity(Vector(X[i][j], Y[i][j], 0))
            u_src[i][j] = V_src.x
            v_src[i][j] = V_src.y
            
            
            V_dblt = doublet_panel.induced_velocity(Vector(X[i][j], Y[i][j], 0))
            u_dblt[i][j] = V_dblt.x
            v_dblt[i][j] = V_dblt.y
    
       
    ax1_src, _ = panel.plot(PanelFixedFrame=False)
    levels = np.linspace(phi_src.min(), phi_src.max(), 20)
    CS=ax1_src.contour(X, Y, phi_src, levels=levels)
    ax1_src.clabel(CS, levels, inline = 1, fmt ='% 1.1f', fontsize = 8)
    
    ax2_src, _  = panel.plot(PanelFixedFrame=False)  
    ax2_src.quiver(X, Y, u_src, v_src, color='m')
        
    ax2_src, _ = panel.plot(PanelFixedFrame=False)
    ax2_src.streamplot(X.T, Y.T, u_src.T, v_src.T)
    
    ax1_dblt, _ = panel.plot(PanelFixedFrame=False)
    levels = np.linspace(phi_dblt.min(), phi_dblt.max(), 20)
    CS=ax1_dblt.contour(X, Y, phi_dblt, levels=levels)
    ax1_dblt.clabel(CS, levels, inline = 1, fmt ='% 1.1f', fontsize = 8)
    
    ax2_dblt, _  = panel.plot(PanelFixedFrame=False)  
    ax2_dblt.quiver(X, Y, u_dblt, v_dblt, color='m')
    x, y = [panel.r[0].x, panel.r[1].x], [panel.r[0].y, panel.r[1].y]
        
    ax2_dblt, _ = panel.plot(PanelFixedFrame=False)
    ax2_dblt.streamplot(X.T, Y.T, u_dblt.T, v_dblt.T)
    
    plt.show()
