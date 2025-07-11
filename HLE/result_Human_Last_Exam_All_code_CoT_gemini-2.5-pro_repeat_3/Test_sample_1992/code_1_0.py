import sympy

def find_qtfp():
    """
    This function solves for the number of Quantum Temporal Fixed Points (QTFPs)
    by setting up and solving the derived algebraic equation.
    """
    # Define the variable 'p' representing the proposition's truth value.
    p = sympy.Symbol('p')

    # The problem of finding QTFPs simplifies to solving a quadratic equation.
    # The equation is derived from equating the forward and backward time-flow operations:
    # sqrt(p**2 + (1-p)**2) = sqrt(p*(1-p) + (1-p)*p)
    # Squaring both sides and simplifying leads to: 4*p**2 - 4*p + 1 = 0.

    # Define the coefficients of the final quadratic equation.
    a = 4
    b = -4
    c = 1

    # Construct the equation object using sympy.
    equation = sympy.Eq(a * p**2 + b * p + c, 0)

    # As requested, we will print the final equation with each number.
    print("The simplified equation for finding Quantum Temporal Fixed Points is:")
    print(f"{a} * p**2 + ({b}) * p + {c} = 0")

    # Solve the equation for p.
    solutions = sympy.solve(equation, p)

    # A proposition's truth value 'p' must be a real number between 0 and 1.
    # We filter the solutions to find the ones that are valid.
    valid_solutions = []
    for s in solutions:
        # Check if the solution is a real number and within the valid probability range [0, 1].
        if s.is_real and 0 <= s <= 1:
            valid_solutions.append(s)

    # The number of valid solutions corresponds to the number of QTFPs.
    num_fixed_points = len(valid_solutions)

    print(f"\nThe unique solution for p is: {valid_solutions[0] if valid_solutions else 'None'}")
    print(f"The number of valid Quantum Temporal Fixed Points is: {num_fixed_points}")

if __name__ == "__main__":
    find_qtfp()