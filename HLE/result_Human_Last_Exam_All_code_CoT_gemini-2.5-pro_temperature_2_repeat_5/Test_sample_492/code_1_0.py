import sympy

def solve_for_c():
    """
    This function calculates the exact value of c, the time of emergence of the
    giant component in the specified random graph model.
    It follows the plan outlined above using symbolic mathematics.
    """
    # Define symbols for time t and arrival time s.
    # We are looking for the critical time c, which we represent with t.
    t, s = sympy.symbols('t s', real=True, positive=True)

    # Step 1: Define the expected degree d(s) of a vertex arriving at time s.
    # In the n->infinity limit, this is found to be d(s) = (t^2 - s^2) / 2.
    d_s = (t**2 - s**2) / 2
    print("Step 1: The expected degree of a vertex arriving at time s in a graph observed at time t is:")
    print(f"d(s) = {d_s}")
    print("-" * 30)

    # The arrival time s of a random vertex is uniformly distributed in [0, t].
    # Its probability density function (pdf) is 1/t.
    pdf_s = 1/t

    # Step 2: Calculate E[D], the average degree of a random vertex.
    # This is the integral of d(s) times the pdf of s, from 0 to t.
    E_D = sympy.integrate(d_s * pdf_s, (s, 0, t))
    print("Step 2: Calculate the average degree E[D].")
    print(f"E[D] = Integral(d(s) * pdf(s), (s, 0, t))")
    print(f"E[D] = Integral(({d_s}) * ({pdf_s}), (s, 0, t)) = {E_D}")
    print("-" * 30)

    # Step 3: Calculate E[D^2], the second moment of the degree distribution.
    # This is the integral of d(s)^2 times the pdf of s, from 0 to t.
    E_D2 = sympy.integrate(d_s**2 * pdf_s, (s, 0, t))
    print("Step 3: Calculate the second moment of the degree distribution E[D^2].")
    print(f"E[D^2] = Integral(d(s)^2 * pdf(s), (s, 0, t))")
    print(f"E[D^2] = Integral(({d_s})**2 * ({pdf_s}), (s, 0, t)) = {E_D2}")
    print("-" * 30)

    # Step 4: Apply the Molloy-Reed criterion for the emergence of a giant component:
    # E[D^2] / E[D] = 2.
    print("Step 4: Set up the equation for the critical time using the Molloy-Reed criterion.")
    # The numbers in the final equation are the coefficients of the symbolic expressions.
    print(f"Criterion: E[D^2] / E[D] = 2")
    print(f"Equation: ({E_D2}) / ({E_D}) = 2")
    
    # Simplify the equation to solve for t
    # (2*t**4/15) / (t**2/3) = (2*t**4 * 3) / (15 * t**2) = 6*t**2 / 15 = 2*t**2 / 5
    simplified_equation = sympy.Eq(sympy.simplify(E_D2 / E_D), 2)
    print(f"Simplified equation: {simplified_equation}")
    print("-" * 30)

    # Step 5: Solve the equation for t to find the constant c.
    solutions = sympy.solve(simplified_equation, t)
    # We are looking for a positive time t.
    c = solutions[0]
    
    print("Step 5: Solve for t to find the critical time c.")
    print(f"The equation {simplified_equation} rearranges to t^2 = 5.")
    print(f"The positive solution for the critical time c is: {c}")

if __name__ == '__main__':
    solve_for_c()
