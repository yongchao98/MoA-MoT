import sympy

def solve_exponent_problem():
    """
    This function solves for the critical exponent p.

    The problem of finding the largest p for which the function I is not in L^p(R^9)
    reduces to analyzing the convergence of an integral. The convergence is determined
    by the slowest decay rate of the function I(a), which occurs on the subspace of
    cubic coefficients.

    The integral over this 4D subspace behaves like:
    ∫ r^(3 - 2p/3) dr, where r is the magnitude of the cubic coefficients.

    This integral diverges if the exponent is >= -1. We are looking for the largest p
    for which this happens, which corresponds to the boundary case.
    """
    
    # We need to solve the equation: 3 - 2p/3 = -1
    # Define the variable and the equation
    p = sympy.Symbol('p')
    
    # These are the numbers from the equation's derivation
    exponent_from_volume = 3
    decay_exponent_numerator = 2
    decay_exponent_denominator = 3
    critical_exponent_value = -1

    # Formulate the equation to be solved
    equation = sympy.Eq(exponent_from_volume - (decay_exponent_numerator * p) / decay_exponent_denominator, critical_exponent_value)
    
    # Solve the equation
    solution = sympy.solve(equation, p)
    p_value = solution[0]

    print("The critical value for the exponent p is determined by the equation:")
    # Print the equation with its constituent numbers
    print(f"{exponent_from_volume} - ({decay_exponent_numerator} * p) / {decay_exponent_denominator} = {critical_exponent_value}")
    
    print("\nSolving for p:")
    step1 = f"{exponent_from_volume} - {critical_exponent_value} = ({decay_exponent_numerator} * p) / {decay_exponent_denominator}"
    print(step1)
    step2 = f"{exponent_from_volume - critical_exponent_value} = ({decay_exponent_numerator} * p) / {decay_exponent_denominator}"
    print(step2)
    step3 = f"{(exponent_from_volume - critical_exponent_value) * decay_exponent_denominator} = {decay_exponent_numerator} * p"
    print(step3)
    step4 = f"p = {((exponent_from_volume - critical_exponent_value) * decay_exponent_denominator)} / {decay_exponent_numerator}"
    print(step4)

    print(f"\nThe solution to the equation is p = {p_value}.")
    
    print(f"\nFor p <= {p_value}, the integral diverges. Therefore, the largest value of p such that the function I is not in L^p(R^9) is {p_value}.")

solve_exponent_problem()
