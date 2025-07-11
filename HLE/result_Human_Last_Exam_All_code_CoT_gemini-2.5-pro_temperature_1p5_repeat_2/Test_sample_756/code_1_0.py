import sys

def solve():
    """
    Finds the maximum value of |b| + |c| for a quadratic ax^2+bx+c
    given that |ax^2+bx+c| <= 1 for all x in [-1, 1].
    """
    print("Let f(x) = ax^2 + bx + c.")
    print("We are given that |f(x)| <= 1 for all x in the interval [-1, 1].")
    print("We want to find the maximum value of |b| + |c|.\n")

    print("Step 1: Express b and c in terms of values of f(x).")
    print("c = f(0)")
    print("b = (f(1) - f(-1)) / 2")
    print("So, |b| + |c| = |(f(1) - f(-1)) / 2| + |f(0)|.\n")

    print("Step 2: Find the conditions for the maximum value.")
    print("To maximize |b| + |c|, we should intuitively choose f(x) to reach the boundaries of the constraint, |f(x)| <= 1.")
    print("Let's choose f(1) and f(-1) to maximize the term |f(1) - f(-1)|. This happens when f(1) and f(-1) have opposite signs and magnitude 1.")
    print("Let's assume f(1) = 1 and f(-1) = -1. This choice maximizes |b|.")
    
    b_val = (1 - (-1)) / 2
    print(f"This gives b = (1 - (-1)) / 2 = {b_val}. So, |b| = {abs(b_val)}.")
    
    print("We now have a polynomial of the form f(x) = ax^2 + x + c.\n")

    print("Step 3: Determine the coefficient 'a' and the polynomial form.")
    print("From f(1) = a + b + c = a + 1 + c = 1, we get a = -c.")
    print("So the polynomial is f(x) = -cx^2 + x + c.")
    print("We also have the constraint |f(0)| = |c| <= 1.\n")

    print("Step 4: Check the condition |f(x)| <= 1 for all x in [-1, 1].")
    print("The extremum (vertex) of this parabola is at x_v = -1 / (2 * (-c)) = 1 / (2c).")
    print("If |c| < 1/2, the vertex is outside [-1, 1], and since f(1)=1 and f(-1)=-1, the condition |f(x)|<=1 is met.")
    print("If |c| >= 1/2, the vertex is inside [-1, 1]. The value at the vertex must also be at most 1 in magnitude.")
    print("The value is f(x_v) = 1/(4c) + c.")
    print("The constraint |f(x_v)| <= 1 implies |1/(4c) + c| <= 1.")
    print("Assuming c > 0, this means 1/(4c) + c <= 1, which simplifies to (2c - 1)^2 <= 0.")
    print("Since a square cannot be negative, the only solution is c = 1/2.\n")
    
    print("Step 5: Conclude the maximum value.")
    print("The analysis shows that for |b|=1, the maximum valid value for |c| is 1/2.")
    print("This gives the maximum value for |b| + |c|.")
    
    b_final = 1.0
    c_final = 0.5
    max_value = abs(b_final) + abs(c_final)

    print("\nThe polynomial that achieves this is f(x) = -0.5*x^2 + x + 0.5.")
    print("For this polynomial: a = -0.5, b = 1, c = 0.5.")
    print("\nThe maximum value is calculated as:")
    print(f"|b| + |c| = |{b_final}| + |{c_final}| = {max_value}")

solve()
<<<1.5>>>