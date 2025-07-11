import sympy

def calculate_final_amplitude():
    """
    This function derives and prints the symbolic expression for the electric field
    amplitude of an electromagnetic wave after passing through a time-varying slab.

    The derivation is based on the following physical principles:
    1. Impedance Matching: The slab's impedance Z(t) = sqrt(mu(t)/epsilon(t)) is
       equal to the vacuum impedance Z0. This means there are no reflections at the
       boundaries.
    2. Frequency Shift: The wave's frequency changes as it propagates through the
       time-varying medium. The frequency at the exit (x=L) is omega_L = omega_0 * exp(-alpha*L/c).
    3. Conservation of Action: For a non-dissipative, time-varying medium, the
       action flux (S/omega) is conserved. Since the Poynting vector S is
       proportional to the amplitude squared (E_amp^2), we have E_amp^2 / omega = constant.

    Combining these principles yields the final amplitude.
    """

    # Define the symbolic variables for the equation
    A = sympy.Symbol('A')          # Initial amplitude of the electric field
    L = sympy.Symbol('L')          # Length of the slab
    alpha = sympy.Symbol('alpha')  # Time-variation parameter of the medium
    c = sympy.Symbol('c')          # Speed of light in vacuum

    # From E_L^2 / omega_L = A^2 / omega_0 and omega_L / omega_0 = exp(-alpha*L/c),
    # we get E_L^2 = A^2 * exp(-alpha*L/c).
    # Taking the square root gives the final amplitude E_L.
    # The numbers in the final equation are -1 and 2 (in the denominator of the exponent).
    
    # Construct the final expression as a string to ensure clear output
    # of all symbols and numbers.
    final_amplitude_expression = f"A * exp(-(alpha * L) / (2 * c))"

    print("The calculated amplitude of the electric field at the rightmost boundary of the slab (E_L) is:")
    print(f"E_L = {final_amplitude_expression}")
    
    print("\nWhere:")
    print(f"  A = Initial amplitude")
    print(f"  L = Length of the slab")
    print(f"  alpha = Time-variation parameter")
    print(f"  c = Speed of light in vacuum")
    print(f"  The number 2 in the denominator of the exponent arises from taking the square root of the power ratio.")

calculate_final_amplitude()