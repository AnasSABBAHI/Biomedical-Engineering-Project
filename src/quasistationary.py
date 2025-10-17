"""
Quasi-stationary regime oxygen diffusion solver.
Solves time-dependent diffusion with periodic boundary conditions.
"""

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from .constants import PhysicalConstants

def solve_quasistationary_diffusion(N, M, L, time, C_a=None, C_b=None, 
                                  C_1=None, omega=None, lambda_param=None):
    """
    Solve quasi-stationary diffusion with time-dependent Dirichlet boundary.
    
    Parameters
    ----------
    N, M : int
        Grid dimensions
    L : float
        Domain length (m)
    time : float
        Current time (s)
    C_a, C_b, C_1 : float
        Concentration parameters (mol/m³)
    omega : float
        Breathing angular frequency (rad/s)
    lambda_param : float
        Screening length (m)
    
    Returns
    -------
    concentration : ndarray
        2D concentration field at given time
    """
    # Default values
    if C_a is None:
        C_a = PhysicalConstants.C_AIR
    if C_b is None:
        C_b = PhysicalConstants.C_BLOOD
    if C_1 is None:
        C_1 = PhysicalConstants.C_REST
    if omega is None:
        omega = PhysicalConstants.OMEGA_REST
    if lambda_param is None:
        lambda_param = PhysicalConstants.LAMBDA_TYPICAL
    
    dx = L / N
    total_points = N * M
    
    # Time-dependent boundary condition
    C_top = C_a - C_b + C_1 * (np.cos(omega * time) - 1)
    
    # Build matrix (same structure as stationary case)
    A = lil_matrix((total_points, total_points))
    B = np.zeros(total_points)
    
    def index(i, j):
        return j * N + i
    
    for i in range(N):
        for j in range(M):
            k = index(i, j)
            
            if 0 < i < N-1 and 0 < j < M-1:
                A[k, index(i+1, j)] = 1
                A[k, index(i-1, j)] = 1
                A[k, index(i, j+1)] = 1
                A[k, index(i, j-1)] = 1
                A[k, k] = -4
                B[k] = 0
                
            elif j == M-1:  # Top boundary (time-dependent Dirichlet)
                A[k, k] = 1
                B[k] = C_top
                
            elif j == 0:  # Bottom boundary (Robin)
                A[k, k] = 1 + dx/lambda_param
                A[k, index(i, 1)] = -1
                B[k] = 0
                
            elif i == 0 and 0 < j < M-1:  # Left boundary (Neumann)
                A[k, k] = 1
                A[k, index(1, j)] = -1
                B[k] = 0
                
            elif i == N-1 and 0 < j < M-1:  # Right boundary (Neumann)
                A[k, k] = 1
                A[k, index(N-2, j)] = -1
                B[k] = 0
    
    A_csr = A.tocsr()
    solution = spsolve(A_csr, B)
    concentration = solution.reshape((M, N)) + C_b
    
    return concentration, C_top

def animate_solution(N=50, M=50, L=0.01, duration=10, fps=10):
    """
    Create animation of quasi-stationary solution.
    
    Parameters
    ----------
    N, M : int
        Grid dimensions
    L : float
        Domain length (m)
    duration : float
        Animation duration (s)
    fps : int
        Frames per second
    
    Returns
    -------
    animation : FuncAnimation
        Matplotlib animation object
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Initial solution
    time = 0
    C, _ = solve_quasistationary_diffusion(N, M, L, time)
    
    x = np.linspace(0, L, N)
    y = np.linspace(0, L, M)
    X, Y = np.meshgrid(x, y)
    
    # Plot 1: Concentration field
    im = ax1.pcolormesh(X, Y, C, shading='auto', cmap='viridis')
    plt.colorbar(im, ax=ax1, label='Concentration (mol/m³)')
    ax1.set_title('Oxygen Concentration Field')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    
    # Plot 2: Time series of flux
    times = []
    fluxes = []
    line, = ax2.plot([], [], 'b-', linewidth=2)
    ax2.set_xlim(0, duration)
    ax2.set_ylim(-2, 2)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Oxygen Flux (mol/s)')
    ax2.set_title('Oxygen Flux vs Time')
    ax2.grid(True)
    
    def update(frame):
        nonlocal times, fluxes
        time = frame / fps
        C, C_top = solve_quasistationary_diffusion(N, M, L, time)
        
        # Update concentration plot
        im.set_array(C.ravel())
        
        # Calculate and update flux
        dx = L / N
        flux = np.sum(C[0, :] / PhysicalConstants.LAMBDA_TYPICAL) * dx
        
        times.append(time)
        fluxes.append(flux)
        
        # Keep only recent data
        if len(times) > 100:
            times = times[-100:]
            fluxes = fluxes[-100:]
        
        line.set_data(times, fluxes)
        
        return im, line
    
    frames = int(duration * fps)
    animation = FuncAnimation(fig, update, frames=frames, 
                            interval=1000/fps, blit=True)
    
    plt.tight_layout()
    return animation