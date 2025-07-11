from sympy import solve, symbols, sqrt

def is_rational_square(num_expr):
    """
    Checks if a symbolic number expression from SymPy is a square of a rational number.
    """
    # The square root of a rational square number must be rational.
    s_root = sqrt(num_expr)
    # SymPy's is_rational is True for rational numbers (integers and fractions).
    return s_root.is_rational

print("We want to determine if 36 + 24*sqrt(2) is a square in the field Q(sqrt(2)).")
print("This is equivalent to finding rational numbers a, b such that (a + b*sqrt(2))^2 = 36 + 24*sqrt(2).")
print("Expanding the left side gives a^2 + 2*b^2 + 2*a*b*sqrt(2).")
print("Comparing coefficients, we get the system of equations:")
print("1) a^2 + 2*b^2 = 36")
print("2) 2*a*b = 24  => b = 12/a")
print("Substituting (2) into (1) gives the equation: a^2 + 2*(144/a^2) = 36")
print("This simplifies to the quartic equation: a^4 - 36*a^2 + 288 = 0.")
print("Let's solve for a^2 by substituting x = a^2:")
print("x^2 - 36*x + 288 = 0\n")

x = symbols('x')
solutions_for_x = solve(x**2 - 36*x + 288, x)

print(f"The solutions for x (which represents a^2) are: {solutions_for_x[0]} and {solutions_for_x[1]}.")

sol1_is_sq = is_rational_square(solutions_for_x[0])
sol2_is_sq = is_rational_square(solutions_for_x[1])

print(f"\nFor 'a' to be rational, a^2 must be the square of a rational number.")
print(f"Is {solutions_for_x[0]} the square of a rational number? {sol1_is_sq}")
print(f"Is {solutions_for_x[1]} the square of a rational number? {sol2_is_sq}")

if not sol1_is_sq and not sol2_is_sq:
    print("\nConclusion: Since a^2 is not the square of a rational number, 'a' cannot be rational.")
    print("This proves that 36 + 24*sqrt(2) is not a square in Q(sqrt(2)), which in turn implies that beta is not a square in K.")
