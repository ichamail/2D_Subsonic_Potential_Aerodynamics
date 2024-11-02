import numpy as np
from scipy.stats import beta
from matplotlib import pyplot as plt


def cosspace(start, end, num_points):
    mean = (start+end)/2
    amp = (end-start)/2
    return mean + amp * np.cos(np.linspace(np.pi, 0, num_points))

def DenserAtBoundaries(start, end, num_points, alpha):
    '''
    Beta distribution
    
    Cumulative distribution function of beta distribution
    
    alpha exists in (-oo, +oo)
    when alpha is 1 => evenly spaced
    when alpha < 1 denser at boundaries
    when alpha > 1 denser at midpoint
    '''
    x = np.linspace(0, 1, num_points)
    a = b = 2-alpha
    return start + beta.cdf(x, a, b) * (end-start)

def DenserAtLeadingEdge(start, end, num_points, factor=1.5):
    """
    Cumulative distribution function of beta distribution
    
    factor > 1
    if factor = 1 evenly spaced
    if factor > 1 denser at leading edge
    (leading edge at end, trailing adge at start)
    """

    x = np.linspace(0, 1, num_points)
    a = 1
    b = factor
    return start + beta.cdf(x, a, b) * (end-start)

def DenserAtTrailingEdge(start, end, num_points, factor=1.5):
    """
    Cumulative distribution function of beta distribution
    
    factor > 1
    if factor = 1 evenly spaced
    if factor > 1 denser at trailing edge
    (leading edge at end, trailing adge at start)
    """

    x = np.linspace(0, 1, num_points)
    a = factor
    b = 1
    return start + beta.cdf(x, a, b) * (end-start)
