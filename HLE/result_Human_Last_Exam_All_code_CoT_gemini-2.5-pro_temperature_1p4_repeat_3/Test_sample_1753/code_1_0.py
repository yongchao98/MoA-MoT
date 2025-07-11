import sympy
from sympy import Symbol, cos, sin, pi, sqrt, integrate, Eq, solve, Rational

def solve_arc_length_problem():
    """
    Solves for the constant 'a' given the parametric equation of an astroid arc
    and its length.
    """
    # Step 1 & 2: Define symbols, equations, and derivatives
    t = Symbol('t')
    x = cos(t)**3
    y = sin(t)**3

    dx_dt = sympy.diff(x, t)
    dy_dt = sympy.diff(y, t)

    # Step 3: Simplify the integrand for arc length
    # (dx/dt)^2 + (dy/dt)^2 simplifies to 9*(sin(t)*cos(t))^2
    # The square root is |3*sin(t)*cos(t)|
    # We will derive the final relation between length L and 'a'
    # L(a) = 3 * a^(2/3)

    # Step 4, 5, 6: Explain the derivation and the relation between L and a.
    print("The arc is defined by the condition 0 <= x <= a.")
    print("This corresponds to two symmetric pieces of the astroid in the first and fourth quadrants.")
    print("The total length L of these two pieces can be expressed as a function of 'a':")
    
    a = Symbol('a')
    length_expr = 3 * a**Rational(2, 3)
    print(f"L(a) = {length_expr}")
    
    # Step 7: Solve for 'a' using the given length
    given_length = Rational(3, 2)
    
    print("\nWe are given that the length of the arc is 3/2.")
    print("We can set up the following equation to solve for 'a':")
    
    # Final Equation Printout
    final_equation = Eq(length_expr, given_length)
    print(f"Equation:  {final_equation.lhs} = {final_equation.rhs}")
    
    # Outputting the numbers in the final equation
    # The equation is 3 * a^(2/3) = 3/2
    num1, num2, num3, num4, num5 = 3, 2, 3, 3, 2
    print(f"The numbers in '{num1} * a^({num2}/{num3}) = {num4}/{num5}' will be used for solving.")

    # Solve the equation for 'a'
    # a^(2/3) = (3/2) / 3 = 1/2
    # a = (1/2)^(3/2)
    solution = solve(final_equation, a)
    
    # The solution will have multiple roots in the complex plane, we need the positive real one.
    positive_real_solution = None
    for sol in solution:
        if sol.is_real and sol > 0:
            positive_real_solution = sol
            break
            
    print("\nSolving for 'a':")
    print(f"a^(2/3) = ({given_length}) / 3 = {given_length/3}")
    print(f"a = ({given_length/3})^(3/2)")

    print("\n----------------------------------------------------")
    print(f"The calculated value of a is: {positive_real_solution}")
    print(f"As a decimal, a is approximately: {positive_real_solution.evalf()}")
    print("----------------------------------------------------")
    
    return positive_real_solution

if __name__ == '__main__':
    final_a = solve_arc_length_problem()
    # The line below is for the final answer block.
    # To represent sqrt(2)/4 in a way that can be parsed.
    # print(f"<<<{final_a}>>>")


solve_arc_length_problem()