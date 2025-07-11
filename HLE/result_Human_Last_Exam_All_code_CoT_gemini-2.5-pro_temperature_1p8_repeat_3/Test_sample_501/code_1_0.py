import sympy as sp

def solve_polymer_force():
    """
    This function derives and prints the force law for a thermally isolated polymer chain.
    
    The force law f(x) gives the restoring force exerted by the polymer when its ends
    are separated by a distance x.

    The formula involves:
    - E0: The kinetic energy of the polymer at zero extension (x=0).
    - x: The separation distance between the polymer ends.
    - n: The number of segments in the polymer chain.
    - l: The length of each segment.
    """

    # Define the symbols used in the equation
    E0 = sp.Symbol('E(0)')
    x, n, l = sp.symbols('x n l', real=True, positive=True)

    # Derived relationship for energy E(x) as a function of extension x for an
    # adiabatic (constant entropy) process.
    exponent = (3 * x**2) / (2 * n**2 * l**2)
    E_x = E0 * sp.exp(exponent)
    
    # The force is the negative derivative of the energy with respect to extension
    # at constant entropy.
    force = -sp.diff(E_x, x)
    
    # Simplify the resulting expression
    force_simplified = sp.simplify(force)

    # --- Construct the output string manually for clarity ---
    # The result from sympy is correct, but we can format it nicely.
    # The force f is a restoring force, hence the negative sign.
    # f = - (pre-factor) * (exponential term)
    # pre-factor = (3 * E(0) * x) / (n^2 * l^2)
    
    print("The derived force law f(x) for the thermally isolated polymer is:")
    
    equation_str = "f(x) = - (3 * E(0) * x / (n**2 * l**2)) * exp(3 * x**2 / (2 * n**2 * l**2))"
    
    # Print the equation with numerical coefficients explicitly shown
    print("\nIn a more readable format:")
    print(f"f(x) = - ( (3 * E(0) * x) / (n^2 * l^2) ) * exp( (3 * x^2) / (2 * n^2 * l^2) )")

    # To be extremely explicit about the numbers as requested:
    numerator_coeff_1 = 3
    numerator_coeff_2 = 3
    denominator_coeff = 2
    
    print(f"\nBreaking down the expression with explicit numbers:")
    print(f"Force = - ( ({numerator_coeff_1} * E(0) * x) / (n^2 * l^2) ) * exp( ({numerator_coeff_2} * x^2) / ({denominator_coeff} * n^2 * l^2) )")
    
    # Return the sympy object for verification if needed, although it's not printed
    return force_simplified

# Execute the function to print the result
solve_polymer_force()

# Final answer in the required format
# Let's write the formula in a simple string for the final answer block.
final_answer_str = "f(x) = -(3*E(0)*x / (n**2 * l**2)) * exp(3*x**2 / (2*n**2 * l**2))"
print(f"\n<<<final_answer>>>\n{final_answer_str}\n<<<final_answer>>>")
