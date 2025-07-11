import math

def calculate_time_dilation():
    """
    Calculates the gravitational time dilation factor 'f' for the Pioneer probe
    and the memory usage 'z' for an efficient C program on the Wuxing architecture.
    """

    # Step 1: Define physical constants
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_SUN = 1.989e30   # Mass of the Sun in kg
    C = 2.99792458e8 # Speed of light in m/s

    # Step 2: Define problem-specific parameters
    M_PANDORA = 2 * M_SUN  # Mass of Pandora is 2 solar masses
    
    # Calculate distance 'r' from the center of Pandora in meters
    R_PANDORA_KM = 20
    D_SURFACE_KM = 60
    r_meters = (R_PANDORA_KM + D_SURFACE_KM) * 1000

    # Step 3: Calculate the Schwarzschild Radius (Rs)
    Rs = (2 * G * M_PANDORA) / (C**2)

    # Step 4: Calculate the time dilation factor 'f'
    # f = sqrt(1 - Rs / r)
    dilation_factor = math.sqrt(1 - (Rs / r_meters))
    
    # Round 'f' to 0.0001 (4 decimal places)
    f_rounded = round(dilation_factor, 4)

    # Step 5: Define the memory usage 'z' in decimal digits (D)
    # The most memory-efficient C program would use 3 variables:
    # - f (frac type): 6D
    # - x (iterator for sqrt, frac type): 6D
    # - i (loop counter, int type): 5D
    z_memory_usage = 6 + 6 + 5

    # Step 6: Print the final equation with numbers, as requested
    print(f"Equation: f = sqrt(1 - {Rs:.0f} / {r_meters:.0f})")
    
    # Step 7: Print the final answer in the required format
    final_answer = f"{f_rounded}:{z_memory_usage}"
    print(final_answer)

calculate_time_dilation()