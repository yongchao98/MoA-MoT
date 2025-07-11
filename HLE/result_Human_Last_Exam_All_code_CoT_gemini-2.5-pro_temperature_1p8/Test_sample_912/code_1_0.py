def solve_magnetic_work():
    """
    This function calculates and prints the symbolic formula for the work done
    by the current source in the described electromechanical cycle.
    """

    # Symbolic representations of the physical quantities
    mu = "μ"
    mu_0 = "μ₀"
    g = "g"
    N = "N"
    w = "w"
    x1 = "x₁"
    x2 = "x₂"
    I1_sq = "I₁²"
    I2_sq = "I₂²"

    # The problem asks for the work done by the current source for each cycle.
    # Based on the derivation, the total work W is given by the formula:
    # W = - (1/2) * (L(x₂) - L(x₁)) * (I₂² - I₁²)
    # where (L(x₂) - L(x₁)) = N² * w/g * (μ - μ₀) * (x₂ - x₁)
    # Substituting this in gives the final expression.

    print("The work done by the current source for each cycle is calculated as:")
    print("W = Final expression")
    print("\nWhere the final expression is constructed as follows:")

    # Build the final equation step-by-step for clarity
    term1 = f"({mu} - {mu_0})"
    term2 = f"(2*{g})"
    term3 = f"{N}²"
    term4 = f"{w}"
    term5 = f"({x2} - {x1})"
    term6 = f"({I2_sq} - {I1_sq})"

    print(f"\nTerm 1 (Permeability difference): {term1}")
    print(f"Term 2 (Gap length factor): {term2}")
    print(f"Term 3 (Turns squared): {term3}")
    print(f"Term 4 (Block width): {term4}")
    print(f"Term 5 (Displacement): {term5}")
    print(f"Term 6 (Current squared difference): {term6}")

    final_equation = f"W = - ( {term1} / {term2} ) * {term3} * {term4} * {term5} * {term6}"

    print("\nCombining these terms gives the final equation:")
    print(final_equation)

solve_magnetic_work()