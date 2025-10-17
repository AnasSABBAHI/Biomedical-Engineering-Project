"""
Boundary condition implementations for diffusion equations.
"""

import numpy as np

class BoundaryCondition:
    """Base class for boundary conditions."""
    
    def __init__(self, value):
        self.value = value
    
    def apply(self, matrix, rhs, indices):
        raise NotImplementedError("Subclasses must implement apply method")

class DirichletBC(BoundaryCondition):
    """Dirichlet boundary condition: u = value."""
    
    def apply(self, matrix, rhs, indices):
        for idx in indices:
            matrix[idx, idx] = 1
            rhs[idx] = self.value

class NeumannBC(BoundaryCondition):
    """Neumann boundary condition: ∂u/∂n = value."""
    
    def apply(self, matrix, rhs, indices, neighbor_indices, dx):
        for idx, neighbor_idx in zip(indices, neighbor_indices):
            matrix[idx, idx] = 1
            matrix[idx, neighbor_idx] = -1
            rhs[idx] = self.value * dx

class RobinBC(BoundaryCondition):
    """Robin boundary condition: ∂u/∂n + alpha * u = value."""
    
    def __init__(self, alpha, value=0):
        self.alpha = alpha
        self.value = value
    
    def apply(self, matrix, rhs, indices, neighbor_indices, dx):
        for idx, neighbor_idx in zip(indices, neighbor_indices):
            matrix[idx, idx] = 1 + self.alpha * dx
            matrix[idx, neighbor_idx] = -1
            rhs[idx] = self.value * dx