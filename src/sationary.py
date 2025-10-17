"""
Stationary regime oxygen diffusion solver.
Solves ΔC = 0 with mixed boundary conditions.
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

from .constants import PhysicalConstants

def solve_stationary_diffusion(N, M, L, C_a=None, C_b=None, lambda_param=None):
    """
    Solve stationary diffusion equation ΔC = 0 with mixed boundary conditions.
    
    Parameters
    ----------
    N, M : int
        Grid dimensions in x and y directions
    L : float
        Domain length in meters
    C_a : float, optional
        Alveolar oxygen concentration (mol/m³). Defaults to PhysicalConstants.C_AIR
    C_b : float, optional  
        Blood oxygen concentration (mol/m³). Defaults to PhysicalConstants.C_BLOOD
    lambda_param : float, optional
        Screening length parameter (m). Defaults to PhysicalConstants.LAMBDA_TYPICAL
    
    Returns
    -------
    concentration : ndarray
        2D array of oxygen concentration (mol/m³)
    """
    # Set default values
    if C_a is None:
        C_a = PhysicalConstants.C_AIR
    if C_b is None:
        C_b = PhysicalConstants.C_BLOOD
    if lambda_param is None:
        lambda_param = PhysicalConstants.LAMBDA_TYPICAL
    
    dx = L / N
    total_points = N * M
    
    # Initialize sparse matrix and RHS vector
    A = lil_matrix((total_points, total_points))
    B = np.zeros(total_points)
    
    def index(i, j):
        """Convert 2D grid indices to 1D array index."""
        return j * N + i
    
    # Build the linear system
    for i in range(N):
        for j in range(M):
            k = index(i, j)
            
            # Interior points: Laplace equation discretization
            if 0 < i < N-1 and 0 < j < M-1:
                A[k, index(i+1, j)] = 1
                A[k, index(i-1, j)] = 1  
                A[k, index(i, j+1)] = 1
                A[k, index(i, j-1)] = 1
                A[k, k] = -4
                B[k] = 0
                
            # Top boundary (Dirichlet): C = C_a - C_b
            elif j == M-1:
                A[k, k] = 1
                B[k] = C_a - C_b
                
            # Bottom boundary (Robin): ∂C/∂n = -C/λ
            elif j == 0:
                A[k, k] = 1 + dx/lambda_param
                A[k, index(i, 1)] = -1
                B[k] = 0
                
            # Left boundary (Neumann): ∂C/∂n = 0
            elif i == 0 and 0 < j < M-1:
                A[k, k] = 1
                A[k, index(1, j)] = -1
                B[k] = 0
                
            # Right boundary (Neumann): ∂C/∂n = 0  
            elif i == N-1 and 0 < j < M-1:
                A[k, k] = 1
                A[k, index(N-2, j)] = -1
                B[k] = 0
                
            # Corner points - average of adjacent boundary conditions
            elif i == 0 and j == 0:  # Bottom-left corner
                A[k, k] = 2 + dx/lambda_param
                A[k, index(1, 0)] = -1
                A[k, index(0, 1)] = -1
                B[k] = 0
            elif i == N-1 and j == 0:  # Bottom-right corner
                A[k, k] = 2 + dx/lambda_param
                A[k, index(N-2, 0)] = -1
                A[k, index(N-1, 1)] = -1
                B[k] = 0
    
    # Solve the linear system
    A_csr = A.tocsr()
    solution = spsolve(A_csr, B)
    
    # Reshape and add blood concentration baseline
    concentration = solution.reshape((M, N)) + C_b
    
    return concentration

def calculate_oxygen_flux(concentration, dx, lambda_param):
    """
    Calculate oxygen flux across boundaries.
    
    Parameters
    ----------
    concentration : ndarray
        2D concentration field
    dx : float
        Grid spacing
    lambda_param : float
        Screening length parameter
    
    Returns
    -------
    flux : float
        Total oxygen flux (mol/s per unit depth)
    """
    # Flux at bottom boundary (Robin condition)
    bottom_flux = np.sum(concentration[0, :] / lambda_param) * dx
    
    return bottom_flux