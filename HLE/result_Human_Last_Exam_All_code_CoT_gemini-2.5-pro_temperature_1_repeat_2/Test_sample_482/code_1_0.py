import sympy

def solve_critical_correlation():
    """
    Uses symbolic mathematics to determine the critical correlation 'C'
    that balances potentiation and depression in the described neural network.
    """
    # Step 1: Define all variables as symbolic objects.
    # C: The covariance (correlation) between corresponding neurons in v and s
    # mu: The average firing rate of input neurons in v and s
    # W_ik_v, W_ik_s: Specific synaptic weights from v_k and s_k to r_i
    # W_i_v, W_i_s: Total synaptic weight from population v and s to neuron r_i
    C, mu = sympy.symbols('C mu')
    W_ik_v, W_ik_s = sympy.symbols('W_ik_v W_ik_s')
    W_i_v, W_i_s = sympy.symbols('W_i_v W_i_s')

    print("--- Derivation of the Critical Correlation ---")
    print("\nThe condition for balancing plasticity from inputs v and s is <r_i * v_k> = <r_i * s_k>.\n")

    # Step 2: Define the expressions for the correlation terms based on the model equations.
    # From the derivation in the thought process, we have:
    # <r_i * v_k> = E[(sum_j W_ij^v v_j + sum_l W_il^s s_l) * v_k]
    # This simplifies based on input statistics (see plan for details).
    corr_r_v = mu * W_ik_v + C * W_ik_s + mu**2 * (W_i_v + W_i_s)

    # By symmetry, swapping v and s gives:
    # <r_i * s_k> = E[(sum_j W_ij^v v_j + sum_l W_il^s s_l) * s_k]
    corr_r_s = C * W_ik_v + mu * W_ik_s + mu**2 * (W_i_v + W_i_s)

    print("Expression for the correlation <r_i * v_k>:")
    print(corr_r_v)
    print("\nExpression for the correlation <r_i * s_k>:")
    print(corr_r_s)

    # Step 3: Set up the balance equation.
    balance_equation = sympy.Eq(corr_r_v, corr_r_s)
    print("\nSetting the two expressions equal gives the balance equation:")
    print(balance_equation)

    # Step 4: Simplify and factor the equation to find the solutions.
    # The mu**2 * (W_i_v + W_i_s) terms on both sides cancel out.
    # We are left with: mu*W_ik_v + C*W_ik_s = C*W_ik_v + mu*W_ik_s
    # Rearranging gives: (C - mu)*W_ik_s - (C - mu)*W_ik_v = 0
    # Factoring gives: (C - mu) * (W_ik_s - W_ik_v) = 0
    
    factored_form = sympy.factor(balance_equation.lhs - balance_equation.rhs)
    print("\nFactoring the simplified equation yields:")
    print(f"{factored_form} = 0")

    # Step 5: Interpret the result.
    print("\nThis equation has two possible solutions:")
    print(f"1. {sympy.Eq(W_ik_v, W_ik_s)}: The weights from both input populations are identical. This corresponds to a non-selective state where the network cannot distinguish between v and s.")
    print(f"2. {sympy.Eq(C - mu, 0)}: This simplifies to C = mu.")

    print("\n--- Conclusion ---")
    print("For the network to develop or maintain selectivity (i.e., for W_ik_v to be different from W_ik_s), the first solution cannot hold.")
    print("Therefore, the only condition that allows for a stable and selective state is the second solution.")
    
    # Final Answer
    final_eq_lhs = C
    final_eq_rhs = mu
    print("\nThe critical amount of correlation 'C' required to balance potentiation and depression is equal to the average input rate 'mu'.")
    print("\nFinal Equation:")
    print(f"{final_eq_lhs} = {final_eq_rhs}")

if __name__ == '__main__':
    solve_critical_correlation()
