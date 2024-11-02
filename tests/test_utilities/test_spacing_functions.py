import numpy as np
from matplotlib import pyplot as plt
from src.utilities import cosspace, DenserAtBoundaries, DenserAtLeadingEdge, DenserAtTrailingEdge

def test_spacing(start, end, num_points, alpha, factor=1.5):


    x = np.linspace(start, end, num_points)
    y1 = np.linspace(start, end, num_points)
    y2 = cosspace(start, end, num_points)
    y3 = DenserAtBoundaries(start, end, num_points, alpha)
    y4 = DenserAtLeadingEdge(start, end, num_points, factor)
    y5 = DenserAtTrailingEdge(start, end, num_points, factor)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    
    ax1.plot(x, y1, "r", label="linear spacing")
    ax1.plot(x, y2, "m", label="cosine spacing")
    ax1.plot(x, y3, "b", label = "prob dens func Beta, α=" + str(alpha))
    ax1.plot(x, y4, "g", label="denser at leading edge")
    ax1.plot(x, y5, "k", label="denser at trailing edge")
    
    ax2.scatter(y1, 2*np.ones_like(y1),s=2, c="r", label="linear spacing")
    ax2.scatter(y2, 1.5*np.ones_like(y2), s=2, c="m", label="cosine spacing")
    ax2.scatter(y3, np.ones_like(y3), s=2, c="b", label = "prob dens func Beta, α=" + str(alpha))
    ax2.scatter(y4, 0.5*np.ones_like(y4), s=2, c="g", label="denser at leading edge")
    ax2.scatter(y5, np.zeros_like(y5), s=2, c="k", label="denser at trailing edge")

    ax1.legend()
    plt.show()
