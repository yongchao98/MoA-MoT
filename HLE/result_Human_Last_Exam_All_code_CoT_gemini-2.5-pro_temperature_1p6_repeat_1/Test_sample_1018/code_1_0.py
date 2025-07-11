import math

def calculate_simplicial_volume_product(g1, g2):
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.
    """
    # The dimension of a surface is 2
    m = 2
    n = 2

    # Step 1: Check if genera are valid for the formula ||Σg|| = 2(2g-2)
    if g1 < 1 or g2 < 1:
        # A surface of genus 0 (a sphere) has simplicial volume 0.
        # This code assumes g >= 1 as per the problem.
        if g1 < 1:
            sv1 = 0
        if g2 < 1:
            sv2 = 0
    
    # Step 2: Calculate simplicial volume for each surface
    # ||Σg|| = 2 * |χ(Σg)| = 2 * (2g - 2) for g >= 1
    sv1 = 2 * (2 * g1 - 2)
    sv2 = 2 * (2 * g2 - 2)

    # Step 3: Calculate the binomial coefficient C(m+n, n)
    # For two surfaces, this is C(2+2, 2) = C(4, 2)
    binomial_coefficient = math.comb(m + n, n)

    # Step 4: Apply the product formula: ||M x N|| = C(m+n, n) * ||M|| * ||N||
    total_sv = binomial_coefficient * sv1 * sv2

    # Print the explanation and final result
    print("The formula for the simplicial volume of an oriented closed surface Σ_g (for g≥1) is ||Σ_g|| = 2(2g - 2).")
    print(f"For g₁ = {g1}, the simplicial volume is ||Σ_{g1}|| = 2(2*{g1} - 2) = {sv1}.")
    print(f"For g₂ = {g2}, the simplicial volume is ||Σ_{g2}|| = 2(2*{g2} - 2) = {sv2}.")
    print("\nThe product formula for simplicial volumes is ||M × N|| = C(m+n, n) * ||M|| * ||N||.")
    print(f"Here, m=2 and n=2, so the binomial coefficient is C(4, 2) = {binomial_coefficient}.")
    print("\nPutting it all together, the final calculation is:")
    print(f"||\u03A3_{g1} \u00D7 \u03A3_{g2}|| = {binomial_coefficient} \u00D7 {sv1} \u00D7 {sv2} = {total_sv}")

# Genera of the surfaces
genus1 = 31
genus2 = 17

calculate_simplicial_volume_product(genus1, genus2)