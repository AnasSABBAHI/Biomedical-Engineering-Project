"""
Physical constants and parameters for oxygen diffusion in pulmonary acinus.
"""

import numpy as np

class PhysicalConstants:
    """Physical constants for oxygen diffusion modeling."""
    
    # Oxygen diffusion coefficient in tissue (m²/s)
    D_O2 = 1.8e-9  
    
    # Oxygen diffusion coefficient in water (m²/s) - for membrane permeability
    D_O2_WATER = 2.1e-9
    
    # Physiological concentrations (mol/m³)
    C_AIR = 8.4          # Oxygen concentration in alveolar air
    C_BLOOD = 5.1e-4     # Oxygen concentration in venous blood
    C_REST = 4.2         # Oxygen concentration in lungs at rest
    
    # Breathing parameters (realistic physiological values)
    BREATHING_RATE_REST = 0.3  # Hz (12-20 breaths/minute at rest)
    BREATHING_RATE_EXERCISE = 0.5  # Hz (during exercise)
    OMEGA_REST = 2 * np.pi * BREATHING_RATE_REST  # rad/s
    OMEGA_EXERCISE = 2 * np.pi * BREATHING_RATE_EXERCISE  # rad/s
    
    # Typical acinus dimensions (meters)
    ACINUS_DIAMETER = 7e-3  # 7 mm
    ACINUS_LENGTH = 10e-3   # 10 mm
    MEMBRANE_THICKNESS = 0.5e-6  # 0.5 μm
    
    # Screening length parameters (m)
    LAMBDA_TYPICAL = 0.28  # 28 cm in meters
    LAMBDA_RANGE = [0.01, 2.0]  # Range for parameter studies
    
    # Membrane permeability parameters
    BETA = 1.0  # Partition coefficient
    TAU = 1.0   # Tortuosity factor