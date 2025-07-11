import math

def calculate_final_amplitude():
    """
    Calculates the amplitude of an EM wave after passing through a time-varying slab.
    
    The final amplitude A_L is given by the formula: A_L = A * exp(- (alpha * L) / (2 * c))
    where:
    A = Initial amplitude of the electric field
    alpha = Time-varying parameter of the slab's material properties
    L = Length of the slab
    c = Speed of light in vacuum
    """
    
    # --- User-defined parameters (example values) ---
    # Initial amplitude of the wave (e.g., in V/m)
    A = 1.0
    # Material parameter alpha (in 1/s)
    alpha = 1.0e8
    # Length of the slab L (in meters)
    L = 1.0
    # Speed of light in vacuum (in m/s)
    c = 299792458.0

    # --- Calculation ---
    exponent = - (alpha * L) / (2 * c)
    final_amplitude = A * math.exp(exponent)

    # --- Output the results ---
    print("--- Calculation of Final Electric Field Amplitude ---")
    print(f"The final amplitude (A_L) is derived from the initial amplitude (A), slab length (L), material parameter (alpha), and the speed of light (c).")
    print("\nDerived Formula:")
    print("A_L = A * exp( - (alpha * L) / (2 * c) )")
    
    print("\nSubstituting the given values:")
    # This part fulfills the requirement to "output each number in the final equation"
    print(f"A_L = {A} * exp( - ({alpha} * {L}) / (2 * {c}) )")
    
    print("\nFinal Result:")
    print(f"The amplitude of the electric field at x=L is: {final_amplitude:.6f} V/m")

if __name__ == "__main__":
    calculate_final_amplitude()