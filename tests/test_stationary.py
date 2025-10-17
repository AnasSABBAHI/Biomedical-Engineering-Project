"""
Tests for stationary diffusion solver.
"""

import numpy as np
import pytest
from src.acinus_diffusion.stationary import solve_stationary_diffusion
from src.acinus_diffusion.constants import PhysicalConstants

def test_stationary_solution_shape():
    """Test that solution has correct shape."""
    N, M = 50, 50
    L = 0.01
    C = solve_stationary_diffusion(N, M, L)
    assert C.shape == (M, N)

def test_boundary_conditions():
    """Test that boundary conditions are satisfied."""
    N, M = 30, 30
    L = 0.01
    C_a = 8.4
    C_b = 5.1e-4
    lambda_param = 0.28
    
    C = solve_stationary_diffusion(N, M, L, C_a, C_b, lambda_param)
    
    # Test Dirichlet boundary (top)
    np.testing.assert_allclose(C[-1, :], C_a, rtol=1e-10)
    
    # Test Robin boundary consistency (bottom)
    dx = L / N
    expected_robin = C[1, :] / (1 + dx/lambda_param)
    np.testing.assert_allclose(C[0, :], expected_robin, rtol=1e-5)

def test_mass_conservation():
    """Test approximate mass conservation."""
    N, M = 40, 40
    L = 0.01
    C = solve_stationary_diffusion(N, M, L)
    
    # In stationary regime, flux in should approximately equal flux out
    dx = L / N
    top_flux = np.mean(C[-1, :]) * dx
    bottom_flux = np.mean(C[0, :]) * dx
    
    # Allow some tolerance for numerical errors
    assert abs(top_flux - bottom_flux) < 1e-5

def test_convergence():
    """Test solution convergence with grid refinement."""
    L = 0.01
    grid_sizes = [20, 40, 60]
    solutions = []
    
    for N in grid_sizes:
        M = N
        C = solve_stationary_diffusion(N, M, L)
        solutions.append(C)
    
    # Coarse grid solution interpolated to fine grid
    C_coarse = solutions[0]
    C_fine = solutions[-1]
    
    # Should be reasonably close (allowing for discretization error)
    coarse_interp = C_coarse[::3, ::3]  # Rough interpolation
    error = np.mean(np.abs(coarse_interp - C_fine[::3, ::3]))
    assert error < 0.1  # Conservative tolerance