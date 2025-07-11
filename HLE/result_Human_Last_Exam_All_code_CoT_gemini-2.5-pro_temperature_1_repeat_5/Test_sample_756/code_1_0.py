from fractions import Fraction

def solve():
    """
    This function demonstrates the solution to the problem.
    It defines the coefficients of the polynomial that gives the maximum value for |b| + |c|,
    calculates this maximum value, and verifies that the polynomial satisfies the given condition.
    """
    
    # The extremal polynomial is f(x) = -(8/9)x^2 + (8/9)x + (7/9).
    # We use the Fraction class for exact arithmetic.
    a = Fraction(-8, 9)
    b = Fraction(8, 9)
    c = Fraction(7, 9)

    print(f"The analysis leads to the extremal polynomial f(x) = ({a})x^2 + ({b})x + ({c})")

    # Calculate |b| + |c| for this polynomial.
    # We can assume b>=0 and c>=0 without loss of generality.
    # The maximum value is b+c.
    max_val = abs(b) + abs(c)

    print("\nFor this polynomial, the value of |b| + |c| is calculated as follows:")
    # The problem asks to output each number in the final equation.
    print(f"|b| + |c| = |{b}| + |{c}| = {abs(b)} + {abs(c)} = {max_val}")
    
    # We must verify that for this polynomial, |f(x)| <= 1 for all x in [-1, 1].
    # Since the parabola opens downwards (a < 0), we check the values at its vertex
    # (which will be the maximum) and at the interval endpoints [-1, 1].
    
    # Define the polynomial function
    def f(x):
        # Ensure x is a Fraction for calculations
        x_frac = Fraction(x)
        return a * x_frac**2 + b * x_frac + c

    # The vertex of the parabola ax^2+bx+c is at x = -b / (2a)
    vertex_x = -b / (2 * a)

    # Calculate the function's value at the vertex and endpoints
    val_at_vertex = f(vertex_x)
    val_at_minus_1 = f(-1)
    val_at_1 = f(1)

    print("\nVerification that |f(x)| <= 1 for x in [-1, 1]:")
    print(f"The vertex is at x = {vertex_x}")
    print(f"The value at the vertex is f({vertex_x}) = {val_at_vertex}. This is the maximum value on the interval.")
    print(f"The value at x = -1 is f(-1) = {val_at_minus_1}. This is the minimum value on the interval.")
    print(f"The value at x = 1 is f(1) = {val_at_1}.")

    if val_at_vertex <= 1 and val_at_minus_1 >= -1:
        print("\nSince the maximum value on [-1, 1] is 1 and the minimum value is -1, the condition |f(x)| <= 1 holds.")
        print(f"Thus, the maximum value of |b| + |c| is indeed {max_val}.")
    else:
        print("\nThe condition |f(x)| <= 1 does not hold for this polynomial.")

solve()