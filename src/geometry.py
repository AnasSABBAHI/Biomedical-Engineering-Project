"""
Geometry creation for different domain types.
"""

import numpy as np

def create_rectangular_domain(N, M, L_x, L_y):
    """
    Create a rectangular computational domain.
    
    Parameters
    ----------
    N, M : int
        Grid dimensions
    L_x, L_y : float
        Domain dimensions in x and y directions (m)
    
    Returns
    -------
    X, Y : ndarray
        2D coordinate arrays
    """
    x = np.linspace(0, L_x, N)
    y = np.linspace(0, L_y, M)
    X, Y = np.meshgrid(x, y)
    return X, Y

def create_deformed_domain(N, M, L_x, L_y, deformation_factor=0.3):
    """
    Create a deformed domain to simulate COPD pathology.
    
    Parameters
    ----------
    N, M : int
        Grid dimensions
    L_x, L_y : float
        Original domain dimensions
    deformation_factor : float
        Amount of deformation (0 = no deformation, 1 = maximum deformation)
    
    Returns
    -------
    X, Y : ndarray
        Deformed coordinate arrays
    mask : ndarray
        Boolean mask indicating valid points
    """
    X, Y = create_rectangular_domain(N, M, L_x, L_y)
    
    # Apply deformation - simulate tissue destruction in COPD
    center_x, center_y = L_x / 2, L_y / 2
    
    # Create elliptical deformation
    rx = L_x * deformation_factor / 2
    ry = L_y * deformation_factor / 2
    
    # Mask for deformed region (simulating destroyed tissue)
    mask = ((X - center_x)**2 / rx**2 + (Y - center_y)**2 / ry**2) <= 1
    
    return X, Y, mask