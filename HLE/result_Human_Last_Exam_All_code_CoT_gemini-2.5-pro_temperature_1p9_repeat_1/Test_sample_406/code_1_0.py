import sys

def get_conditions_for_tripled_fixed_point():
    """
    This script explains and prints the conditions for the existence of an
    FGH-tripled fixed point based on the Banach Fixed-Point Theorem applied
    to a product space.
    """

    # We need to handle unicode characters for mathematical symbols
    if sys.stdout.encoding != 'UTF-8':
        sys.stdout.reconfigure(encoding='UTF-8')

    print("--- Conditions for an FGH-Tripled Fixed Point ---")
    print("\nStep 1: Define the Mathematical Setup")
    print("Let (X, d_X), (Y, d_Y), and (Z, d_Z) be metric spaces.")
    print("Let the functions F, G, and H be defined as follows (correcting the second 'G' to 'H'):")
    print("  F: X × Y × Z → X")
    print("  G: Y × X × Y → Y")
    print("  H: Z × Y × X → Z")

    print("\nStep 2: Define the FGH-Tripled Fixed Point")
    print("A point (x, y, z) ∈ X × Y × Z is an FGH-tripled fixed point if it satisfies the system of equations:")
    print("  x = F(x, y, z)")
    print("  y = G(y, x, y)")
    print("  z = H(z, y, x)")

    print("\nStep 3: Define a Single Mapping on the Product Space")
    print("We can combine these into a single mapping T on the product space M = X × Y × Z.")
    print("Define the metric d_M on M as: d_M((x₁, y₁, z₁), (x₂, y₂, z₂)) = d_X(x₁, x₂) + d_Y(y₁, y₂) + d_Z(z₁, z₂).")
    print("Define the mapping T: M → M as:")
    print("  T(x, y, z) = (F(x, y, z), G(y, x, y), H(z, y, x))")
    print("A fixed point of T, where T(x, y, z) = (x, y, z), is precisely an FGH-tripled fixed point.")

    print("\nStep 4: Apply the Banach Fixed-Point Theorem")
    print("The Banach Fixed-Point Theorem states that if T is a contraction on a complete metric space M, it has a unique fixed point.")
    print("This leads to the following conditions:")

    print("\n--- DERIVED CONDITIONS ---")
    print("\nCONDITION A: Completeness")
    print("The metric spaces (X, d_X), (Y, d_Y), and (Z, d_Z) must be complete.")

    print("\nCONDITION B: Individual Contraction Inequalities")
    print("There must exist non-negative constants (k₁, k₂, k₃, l₁, l₂, l₃, m₁, m₂, m₃) such that for all x₁,x₂ ∈ X, y₁,y₂ ∈ Y, z₁,z₂ ∈ Z:")
    print("  1. d_X(F(x₁, y₁, z₁), F(x₂, y₂, z₂))  ≤  k₁·d_X(x₁, x₂) + k₂·d_Y(y₁, y₂) + k₃·d_Z(z₁, z₂)")
    print("  2. d_Y(G(y₁, x₁, y'₁), G(y₂, x₂, y'₂))  ≤  l₁·d_Y(y₁, y₂) + l₂·d_X(x₁, x₂) + l₃·d_Y(y'₁, y'₂)")
    print("  3. d_Z(H(z₁, y₁, x₁), H(z₂, y₂, x₂))  ≤  m₁·d_Z(z₁, z₂) + m₂·d_Y(y₁, y₂) + m₃·d_X(x₁, x₂)")

    print("\nCONDITION C: Combined Contraction Condition")
    print("To ensure the overall mapping T is a contraction, the sum of coefficients affecting each component metric must be less than 1.")
    print("The sufficient condition is that the maximum of these sums is less than 1.")
    print("\nThe final equation for the condition is:")
    print("  max(k₁ + l₂ + m₃,  k₂ + l₁ + l₃ + m₂,  k₃ + m₁)  <  1")
    print("\nIf all three conditions (A, B, C) are met, a unique FGH-tripled fixed point exists.")


if __name__ == '__main__':
    get_conditions_for_tripled_fixed_point()