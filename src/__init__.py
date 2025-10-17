"""
Acinus Oxygen Diffusion Model

A computational model for simulating oxygen diffusion in pulmonary acinus
using finite difference methods with various boundary conditions.
"""

__version__ = "0.1.0"
__author__ = "Biomedical Engineering Group"

from .stationary import solve_stationary_diffusion
from .quasistationary import solve_quasistationary_diffusion, animate_solution
from .boundary_conditions import DirichletBC, NeumannBC, RobinBC
from .geometry import create_rectangular_domain, create_deformed_domain
from .constants import PhysicalConstants
from .visualization import plot_concentration_field, plot_oxygen_flux

__all__ = [
    'solve_stationary_diffusion',
    'solve_quasistationary_diffusion',
    'animate_solution',
    'DirichletBC',
    'NeumannBC',
    'RobinBC',
    'create_rectangular_domain',
    'create_deformed_domain',
    'PhysicalConstants',
    'plot_concentration_field',
    'plot_oxygen_flux',
]