def solve_poincare_hopf_boundary():
    """
    Calculates the minimum number of zeros of a vector field on a
    compact manifold M with a non-empty boundary ∂M.

    The formula used is min(|χ(M)|, |χ(M) - χ(∂M)|).

    This example uses the 3-ball, where χ(M) = 1 and χ(∂M) = 2.
    """

    # Euler characteristics for the example M = 3-ball
    chi_M = 1
    chi_dM = 2

    print("This script calculates the least number of zeros for a vector field on a compact manifold M with boundary ∂M.")
    print(f"We will use the example of a 3-ball, where χ(M) = {chi_M} and χ(∂M) = {chi_dM}.")
    print("-" * 20)
    print("The formula for the minimum number of zeros is: min(|χ(M)|, |χ(M) - χ(∂M)|)")
    
    # Calculate the two terms in the minimum
    term1 = abs(chi_M)
    term2 = abs(chi_M - chi_dM)
    
    # Calculate the final result
    result = min(term1, term2)

    # Print the equation with the numbers filled in, as requested.
    # The final equation is min(|1|, |1 - 2|) = 1
    print("\nEvaluating the expression:")
    print(f"The first term is |χ(M)| = |{chi_M}| = {term1}")
    print(f"The second term is |χ(M) - χ(∂M)| = |{chi_M} - {chi_dM}| = |{chi_M - chi_dM}| = {term2}")
    print(f"The result is min({term1}, {term2}) = {result}")
    print("-" * 20)
    print(f"The least number of zeros a vector field can have on this manifold is {result}.")

# Run the solver
solve_poincare_hopf_boundary()
