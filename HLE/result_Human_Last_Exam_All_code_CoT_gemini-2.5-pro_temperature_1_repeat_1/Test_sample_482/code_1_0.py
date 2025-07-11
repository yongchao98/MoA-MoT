import sympy as sp

def solve_critical_correlation():
    """
    This function determines the critical amount of correlation between
    input populations v and s that balances potentiation and depression.

    The balance leads to a critical point in the dynamics of weight selectivity.
    Based on the mathematical derivation, the evolution of selectivity (the difference
    between weights W^v and W^s) is governed by a factor (μ - C), where:
    - μ is the average firing rate of the input neurons.
    - C is the correlation between co-located input neurons, defined as <v_k * s_k>.

    The critical point occurs when this factor is zero, marking the transition
    between a regime where selectivity grows and a regime where it decays.
    """

    # Define symbolic variables for the parameters in the equation
    # mu represents the average activation rate of input populations v and s
    mu = sp.Symbol('μ')
    
    # C represents the correlation between the two input populations, <v_k * s_k>
    C = sp.Symbol('C')

    # The critical condition is when the factor governing selectivity dynamics is zero
    critical_equation = sp.Eq(mu - C, 0)

    # Solve the equation for the critical correlation C
    # The result from solve is a list, so we extract the first element
    critical_C_solution = sp.solve(critical_equation, C)[0]

    # Create the final equation object for printing
    final_equation = sp.Eq(C, critical_C_solution)

    # Output the final equation and explain its components
    print("The critical amount of correlation is given by the equation:")
    # We use f-string formatting to present the equation clearly
    print(f"{final_equation.lhs} = {final_equation.rhs}")
    print("\nWhere:")
    print(f"  {final_equation.lhs} is the 'critical amount of correlation'. It is defined as the expectation value of the product of the activities of co-located neurons from the two input populations, i.e., C = <v_k * s_k>.")
    print(f"  {final_equation.rhs} is the average rate of activation of the input neurons, i.e., μ = <v_k> = <s_k>.")

# Run the function to get the answer
solve_critical_correlation()