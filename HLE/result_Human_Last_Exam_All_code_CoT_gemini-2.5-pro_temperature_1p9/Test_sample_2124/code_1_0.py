import sys

def solve_yukawa_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The calculation follows these steps:
    1. Define relationships between the required counter-terms.
    2. Use results from one-loop calculations for the fermion self-energy to find the coefficients
       for delta Z_x and delta Z_m.
    3. Use the known result for the vertex counter-term delta Z_1 in Yukawa theory.
    4. Assemble the final ratio R.
    """

    # Let C be the common divergent factor C = g^2 / (32 * pi^2 * epsilon).
    # We work with the numerical coefficients of C for each counter-term.

    # Step 1: Calculate fermion field and mass counter-terms from the fermion self-energy diagram.
    # From the one-loop fermion self-energy, we find:
    # delta Z_x is the counter-term for the kinetic part (i*gamma*p).
    # delta Z_m is the counter-term for the mass part (M_x).
    
    # delta Z_x = -1 * C
    delta_Zx_factor = -1.0
    # delta Z_m = 3 * C
    delta_Zm_factor = 3.0

    # Step 2: Use the known relation for the vertex counter-term in Yukawa theory.
    # It is a known result that at one-loop, delta Z_1 = 2 * delta Z_x.
    delta_Z1_factor = 2.0 * delta_Zx_factor

    # Step 3: Define the terms appearing in the ratio R.
    # The ratio is R = delta Z_x / (delta Z_g + delta Z_mx)

    # delta Z_g is defined as delta Z_1 - delta Z_x - 0.5 * delta Z_phi
    # The problem states to assume delta Z_phi = 0.
    delta_Z_phi_factor = 0
    delta_Zg_factor = delta_Z1_factor - delta_Zx_factor - 0.5 * delta_Z_phi_factor
    
    # delta Z_mx is defined as delta Z_m + delta Z_x
    delta_Zmx_factor = delta_Zm_factor + delta_Zx_factor

    # Step 4: Compute the ratio R by substituting the factors.
    numerator_factor = delta_Zx_factor
    denominator_factor = delta_Zg_factor + delta_Zmx_factor
    
    # Perform the final calculation
    R = numerator_factor / denominator_factor

    # Print the step-by-step derivation of the final equation
    print("Calculating the ratio R = delta Z_x / (delta Z_g + delta Z_mx).")
    print("The counter-terms are expressed as coefficients of a common factor C = g^2 / (32 * pi^2 * epsilon).")
    print("-" * 50)
    
    print(f"1. Fermion field counter-term: delta Z_x = {delta_Zx_factor} * C")
    
    print("\n2. Yukawa coupling counter-term (with delta Z_phi = 0):")
    print(f"   delta Z_g = delta Z_1 - delta Z_x = (2 * delta Z_x) - delta Z_x = delta Z_x")
    print(f"   delta Z_g = {delta_Zg_factor} * C")

    print("\n3. Fermion mass term counter-term:")
    print(f"   delta Z_mx = delta Z_m + delta Z_x = ({delta_Zm_factor})*C + ({delta_Zx_factor})*C")
    print(f"   delta Z_mx = {delta_Zmx_factor} * C")

    print("-" * 50)
    print("Now, substitute these into the expression for R:")
    
    # The final equation showing each number
    print(f"R = ({delta_Zx_factor} * C) / ( ({delta_Zg_factor} * C) + ({delta_Zmx_factor} * C) )")
    print("The common factor C cancels out:")
    print(f"R = {delta_Zx_factor} / ( {delta_Zg_factor} + {delta_Zmx_factor} )")
    print(f"R = {numerator_factor} / {denominator_factor}")
    print(f"\nThe final result is:")
    print(f"R = {R}")


if __name__ == "__main__":
    solve_yukawa_ratio()
