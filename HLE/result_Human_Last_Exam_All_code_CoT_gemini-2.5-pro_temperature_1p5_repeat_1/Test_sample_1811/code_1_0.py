import sympy

def solve_zeros_on_manifold():
    """
    Calculates and prints the formula for the minimum number of zeros of a vector field
    on a compact manifold M with a non-empty boundary.
    """
    # Define symbolic variables for the Euler characteristics.
    chi_M = sympy.Symbol('χ(M)')
    chi_dM = sympy.Symbol('χ(∂M)')

    # The formula for the least number of zeros.
    # The result depends only on the Euler characteristic of the boundary.
    # We use Abs for the absolute value, as is standard in sympy.
    min_zeros_formula = sympy.Abs(chi_dM / 2)

    # Print the explanation and the final formula.
    print("Let M be a compact manifold with a non-empty boundary ∂M.")
    print("The Euler characteristic of the boundary, χ(∂M), is always an even integer.")
    print("The least number of zeros a vector field can have on M is given by the formula:")
    
    # Using sympy.pretty_print for a nicer output of the formula.
    # The prompt asked to output each number, so we construct the final string manually.
    final_equation_str = f"Minimum Zeros = | {chi_dM} / 2 |"
    print(final_equation_str)

# Execute the function to get the answer.
solve_zeros_on_manifold()

# Final answer expression to be returned
final_answer = "|χ(∂M)/2|"