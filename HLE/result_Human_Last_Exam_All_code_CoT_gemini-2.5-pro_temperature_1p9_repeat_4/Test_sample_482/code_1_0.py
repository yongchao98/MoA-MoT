import sympy

def solve_critical_correlation():
    """
    This function outlines the derivation for the critical correlation 'c'
    and prints the resulting equation.

    The derivation is based on analyzing the fixed point of the weight dynamics.
    The condition for the stable fixed point is that the average weight change for both pathways is zero.
    <dW_v/dt> = 0  and  <dW_s/dt> = 0
    
    This leads to the relation <r_i*v_k> = <r_i*s_l> for any neurons k and l.
    
    Assuming:
    1. Input populations have equal size: N_v = N_s = N
    2. Inputs within a population are uncorrelated.
    3. Correlation between a neuron in v and a neuron in s is uniform and given by the coefficient 'c'.
    4. The system settles into a state with uniform weights w_v and w_s for each pathway.

    Expanding the relation <r_i*v_k> = <r_i*s_l> with the linear neuron model (r_i = w_v * sum(v) + w_s * sum(s))
    and the given input statistics leads to the following equation for the fixed-point weights:

    w_v * (1 - c * N) = w_s * (1 - c * N)

    This equation governs the balance between the two pathways.

    - If (1 - c * N) is not zero, the only solution is w_v = w_s. This is the non-selective state, where
      both pathways have equal strength.

    - If (1 - c * N) = 0, the equation becomes 0 = 0. The constraint w_v = w_s is lifted. The weights
      can evolve independently, allowing for selectivity to develop. This is the "critical point".

    The critical correlation 'c' is the value that satisfies this condition.
    """

    # We need to solve the equation: 1 - c * N = 0 for c.
    # Let's define our symbols
    N = sympy.Symbol('N', positive=True, integer=True)
    c = sympy.Symbol('c')
    
    # The equation to solve
    equation = sympy.Eq(1 - c * N, 0)
    
    # Solve for c
    critical_c_solution = sympy.solve(equation, c)
    
    # The result is a list, so we take the first element
    c_critical = critical_c_solution[0]
    
    # Now, let's print the numbers in the final equation: c = 1 / N
    numerator = 1
    denominator = N
    
    print("Derivation leads to the equation: 1 - c * N = 0")
    print("Solving for 'c' gives the critical amount of correlation.")
    print("\nThe final equation for the critical correlation is C_critical = numerator / denominator.")
    print(f"Where the numerator is: {numerator}")
    print(f"And the denominator is the number of neurons per input layer: {denominator}")
    
    print("\n--- Final Answer ---")
    final_eq = sympy.Eq(sympy.Symbol('C_critical'), c_critical)
    # The pretty print function displays the equation in a formatted way.
    sympy.pprint(final_eq)

solve_critical_correlation()