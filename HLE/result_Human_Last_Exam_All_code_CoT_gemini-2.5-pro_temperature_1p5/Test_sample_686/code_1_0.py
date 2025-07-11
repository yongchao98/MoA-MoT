import sympy

def solve_ising_susceptibility():
    """
    This function derives and prints the symbolic expression for the magnetic susceptibility
    of the Ising model on a sparse random graph as described in the problem.
    """

    # Define symbolic variables
    beta = sympy.Symbol('beta')  # Inverse temperature
    c = sympy.Symbol('c')        # Connectivity
    J = sympy.Symbol('J')        # Homogeneous coupling
    m0 = sympy.Symbol('m0')      # Magnetization at site 0
    mc = sympy.Symbol('mc')      # Cavity magnetization

    # T is the shorthand for tanh(beta*J)
    T = sympy.Symbol('T')

    # M is the propagation factor for the perturbation
    M = sympy.Symbol('M')

    # N is the constant given in the hint
    N = sympy.Symbol('N')

    # The derivation steps outlined in the plan lead to the following expression for chi.
    # Step 1-4: The correlation C_l is found to be C_l = (1 - m0**2) * M**l
    # Step 5-6: Substituting this into the sum for chi gives:
    # chi = beta * Sum(c * (c-1)**(l-1) * (1 - m0**2) * M**l) for l from 1 to oo
    # This is a geometric series which sums to:
    # chi = (beta * c * (1 - m0**2) * M) / (1 - (c - 1) * M)

    # Step 7: Express the result using N = beta * c * (1 - m0**2) / (c - 1)
    # chi = N * (c - 1) * M / (1 - (c - 1) * M)
    # Let's create the final expression.

    final_chi_expr = N * (((c - 1) * M) / (1 - (c - 1) * M))

    # We will print the final equation for chi, including the definitions for the symbols.
    # The user asked to output each number in the final equation. Here, the numbers are 1 and implicitly the powers.

    print("The magnetic susceptibility is given by:")
    print(f"chi = {final_chi_expr}")
    print("\nwhere the terms are defined as:")
    print(f"N = {beta * c * (1 - m0**2) / (c - 1)}")
    print(f"M = {T * (1 - mc**2) / (1 - T**2 * mc**2)}")
    print(f"T = tanh({beta}*{J})")
    print("\nThe magnetizations m0 and mc are determined by the self-consistent equations at zero external field:")
    print(f"mc = tanh(({c} - 1)*atanh({T}*{mc}))")
    print(f"m0 = tanh({c}*atanh({T}*{mc}))")

# Execute the function to print the result
solve_ising_susceptibility()

# Final answer format requested by the user.
# The question is symbolic, so the final answer is the derived formula.
# We present it in the most simplified form based on the derivation.
final_answer_string = "N * ((c - 1) * M) / (1 - (c - 1) * M)"
print(f"\n<<<{final_answer_string}>>>")