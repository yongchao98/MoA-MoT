import math

def display_inner_product_formula():
    """
    This function calculates and displays the formula for the inner product
    (phi, D_phi phi) in finite-temperature field theory for a neutral scalar field.
    """
    
    # Unicode characters for mathematical symbols to make the output clearer
    phi = "\u03D5"
    sum_symbol = "\u2211"
    integral_symbol = "\u222B"
    pi_symbol = "\u03C0"
    omega = "\u03C9"
    
    # --- The numbers and exponents in the equation ---
    num_two_in_denominator = 2
    num_two_in_matsubara = 2
    exponent_field = 2
    exponent_omega = 2
    exponent_k = 2
    exponent_m = 2
    
    # The variable 'd' represents the number of spatial dimensions.
    # The variable 'n' is the integer index for the Matsubara frequencies.

    # --- Constructing the equation string ---
    
    # Left-hand side of the equation
    lhs = f"({phi}, D{phi} {phi})"

    # Right-hand side of the equation
    # T * Sum_{n} Integral[d^d k / (2*pi)^d] |phi_n(k)|^2 * (omega_n^2 + k^2 + m^2)
    rhs = (f"T * {sum_symbol}\u2099  {integral_symbol} [d\u1d48k / ({num_two_in_denominator}{pi_symbol})\u1d48] |{phi}\u2099(k)|\u00B2 "
           f"* ({omega}\u2099\u00B2 + k\u00B2 + m\u00B2)")
    
    # Definition of the Matsubara frequency for a bosonic field
    matsubara_definition = f"where {omega}\u2099 = {num_two_in_matsubara} * {pi_symbol} * n * T"

    # --- Printing the final result ---
    
    print("For a neutral scalar field at finite temperature T, the inner product (phi, D_phi phi)")
    print("is expressed in the Fourier basis (Matsubara frequencies and spatial momentum) as follows:")
    print("-" * 70)
    print(f"{lhs} = {rhs}")
    print("\n" + matsubara_definition)
    print("-" * 70)
    
    # As requested, outputting each number and exponent in the final equation.
    print("\n breakdown of the numerical components and exponents:")
    print(f"  - In the integral's denominator, the number in front of pi ({pi_symbol}) is: {num_two_in_denominator}")
    print(f"  - In the definition of the Matsubara frequency ({omega}\u2099), the integer coefficient is: {num_two_in_matsubara}")
    print(f"  - The exponent of the field amplitude |{phi}\u2099(k)| is: {exponent_field}")
    print(f"  - The exponent of the Matsubara frequency {omega}\u2099 is: {exponent_omega}")
    print(f"  - The exponent of the momentum magnitude k is: {exponent_k}")
    print(f"  - The exponent of the mass m is: {exponent_m}")

if __name__ == "__main__":
    display_inner_product_formula()