"""
Tests for geometry creation functions.
"""

import numpy as np
from src.acinus_diffusion.geometry import create_rectangular_domain, create_deformed_domain

def test_rectangular_domain():
    """Test rectangular domain creation."""
    N, M = 50, 60
    L_x, L_y = 0.01, 0.02
    X, Y = create_rectangular_domain(N, M, L_x, L_y)
    
    assert X.shape == (M, N)
    assert Y.shape == (M, N)
    assert np.allclose(X[0, :], np.linspace(0, L_x, N))
    assert np.allclose(Y[:, 0], np.linspace(0, L_y, M))

def test_deformed_domain():
    """Test deformed domain creation."""
    N, M = 40, 40
    L_x, L_y = 0.01, 0.01
    X, Y, mask = create_deformed_domain(N, M, L_x, L_y, deformation_factor=0.3)
    
    assert X.shape == (M, N)
    assert Y.shape == (M, N)
    assert mask.shape == (M, N)
    assert mask.dtype == bool
    
    # Deformation should create some masked points
    assert np.any(mask) and not np.all(mask)

def test_deformation_factor_range():
    """Test deformation factor boundaries."""
    N, M = 20, 20
    L_x, L_y = 0.01, 0.01
    
    # No deformation
    _, _, mask_zero = create_deformed_domain(N, M, L_x, L_y, 0)
    # Maximum deformation
    _, _, mask_max = create_deformed_domain(N, M, L_x, L_y, 1)
    
    # More deformation should mask more points
    assert np.sum(mask_max) >= np.sum(mask_zero)