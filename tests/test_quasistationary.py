"""
Tests for quasi-stationary diffusion solver.
"""

import numpy as np
import pytest
from src.acinus_diffusion.quasistationary import solve_quasistationary_diffusion

def test_quasistationary_time_dependence():
    """Test that solution changes with time."""
    N, M = 30, 30
    L = 0.01
    times = [0, 1, 2]
    solutions = []
    
    for t in times:
        C, _ = solve_quasistationary_diffusion(N, M, L, t)
        solutions.append(C)
    
    # Solutions at different times should be different
    assert not np.allclose(solutions[0], solutions[1])
    assert not np.allclose(solutions[0], solutions[2])

def test_periodic_behavior():
    """Test periodic behavior matches breathing frequency."""
    N, M = 20, 20
    L = 0.01
    omega = 2 * np.pi * 0.3  # 0.3 Hz
    
    # Solutions at t and t + period should be similar
    t1 = 0
    t2 = 1/0.3  # One period later
    
    C1, _ = solve_quasistationary_diffusion(N, M, L, t1, omega=omega)
    C2, _ = solve_quasistationary_diffusion(N, M, L, t2, omega=omega)
    
    # Should be very close (allowing for numerical errors)
    np.testing.assert_allclose(C1, C2, rtol=1e-2)

def test_boundary_amplitude():
    """Test that boundary condition amplitude is correct."""
    N, M = 25, 25
    L = 0.01
    C_a = 8.4
    C_b = 5.1e-4
    C_1 = 4.2
    omega = 1.0
    
    # Test at specific times
    t_max = 0  # cos(0) = 1 → maximum concentration
    t_min = np.pi / omega  # cos(π) = -1 → minimum concentration
    
    C_max, C_top_max = solve_quasistationary_diffusion(N, M, L, t_max, 
                                                      C_a, C_b, C_1, omega)
    C_min, C_top_min = solve_quasistationary_diffusion(N, M, L, t_min, 
                                                      C_a, C_b, C_1, omega)
    
    expected_max = C_a - C_b + C_1 * (1 - 1)  # = C_a - C_b
    expected_min = C_a - C_b + C_1 * (-1 - 1)  # = C_a - C_b - 2*C_1
    
    np.testing.assert_allclose(C_top_max, expected_max, rtol=1e-10)
    np.testing.assert_allclose(C_top_min, expected_min, rtol=1e-10)