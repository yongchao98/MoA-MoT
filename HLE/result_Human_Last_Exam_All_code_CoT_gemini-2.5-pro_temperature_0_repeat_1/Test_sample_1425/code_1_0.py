import sympy as sp

def solve_partition_function():
    """
    This function derives and prints the partition function Z for a system
    with Hamiltonian H = -μN in the grand canonical ensemble.
    """

    # Define the symbols for our equation
    Z, H_op, N_op, mu, beta = sp.symbols('Z H_op N_op mu beta')
    N = sp.Symbol('N', integer=True, nonnegative=True)
    Tr = sp.Function('Tr')
    exp = sp.exp

    print("Derivation of the Partition Function Z:\n")

    # Step 1: State the general formula
    print("Step 1: The general formula for the grand canonical partition function is:")
    general_formula = sp.Eq(Z, Tr(exp(-beta * (H_op - mu * N_op))))
    print(sp.pretty(general_formula, use_unicode=True))
    print("-" * 40)

    # Step 2: Substitute the given Hamiltonian
    print("Step 2: We substitute the given Hamiltonian H_op = -μ * N_op:")
    H_specific = -mu * N_op
    substituted_formula = general_formula.subs(H_op, H_specific)
    print(sp.pretty(substituted_formula, use_unicode=True))
    print("-" * 40)

    # Step 3: Simplify the exponent
    print("Step 3: Simplify the expression inside the exponent:")
    print("  -β * ((-μ * N_op) - (μ * N_op))")
    print("= -β * (-2 * μ * N_op)")
    print("= 2 * β * μ * N_op")
    simplified_formula = sp.Eq(Z, Tr(exp(2 * beta * mu * N_op)))
    print("\nThis gives the simplified partition function:")
    print(sp.pretty(simplified_formula, use_unicode=True))
    print("-" * 40)

    # Step 4: Evaluate the trace
    print("Step 4: We evaluate the trace (Tr) by summing over the particle number states |N>.")
    print("In this basis, the operator N_op becomes the number N. The trace becomes a sum:")
    sum_expression = sp.Sum(exp(2 * beta * mu * N), (N, 0, sp.oo))
    trace_as_sum = sp.Eq(Z, sum_expression)
    print(sp.pretty(trace_as_sum, use_unicode=True))
    print("-" * 40)

    # Step 5: Sum the geometric series
    print("Step 5: This is a geometric series Σ r^N with r = exp(2*β*μ).")
    print("For the series to converge, we require |r| < 1, which means μ < 0.")
    print("The sum of the series is 1 / (1 - r).")
    final_result = 1 / (1 - exp(2 * beta * mu))
    final_equation = sp.Eq(Z, final_result)

    print("\nThe final result for the partition function Z is:")
    print(sp.pretty(final_equation, use_unicode=True))
    print("\nAs requested, the numbers in the final equation are:")
    print(f"  - The numerator: 1")
    print(f"  - The constant in the denominator: 1")
    print(f"  - The factor in the exponent: 2")

if __name__ == '__main__':
    solve_partition_function()