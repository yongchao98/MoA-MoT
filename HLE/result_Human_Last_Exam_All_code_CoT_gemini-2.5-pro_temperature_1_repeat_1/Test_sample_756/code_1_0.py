import sys

def solve():
    """
    This function solves the problem of finding the maximum value of |b| + |c|
    for a quadratic polynomial f(x) = ax^2 + bx + c, given |f(x)| <= 1 for x in [-1, 1].
    The function prints the logical steps of the solution and the final answer.
    """
    print("The problem is to find the maximum value of |b| + |c| for a quadratic f(x) = ax^2 + bx + c,")
    print("given that |f(x)| <= 1 for all x in the interval [-1, 1].")
    
    print("\nStep 1: Express b and c in terms of function values.")
    print("From f(x) = ax^2 + bx + c:")
    print("c = f(0)")
    print("b = (f(1) - f(-1)) / 2")

    print("\nStep 2: Find bounds on |b| and |c|.")
    print("Using the given condition |f(x)| <= 1:")
    print("|c| = |f(0)| <= 1")
    print("|b| = |f(1) - f(-1)| / 2 <= (|f(1)| + |f(-1)|) / 2 <= (1 + 1) / 2 = 1")
    print("This gives a loose upper bound of |b| + |c| <= 2.")

    print("\nStep 3: Investigate an extremal case.")
    print("The maximum is likely achieved when constraints are tight. Let's assume |b| is at its maximum, |b|=1.")
    print("Let b = 1. This requires f(1) = 1 and f(-1) = -1.")
    
    print("\nStep 4: Characterize the polynomial for b=1.")
    print("If b=1, f(1)=1, and f(-1)=-1, the polynomial must be of the form f(x) = a(x^2 - 1) + x.")
    print("For this polynomial, c = f(0) = a(0 - 1) + 0 = -a.")
    
    print("\nStep 5: Find the constraint on 'a'.")
    print("For the condition |f(x)| <= 1 to hold for f(x) = a(x^2 - 1) + x, it can be shown that we must have |a| <= 0.5.")
    
    print("\nStep 6: Maximize the target expression.")
    print("We want to maximize |b| + |c|. Substituting our findings:")
    print("|b| + |c| = |1| + |-a| = 1 + |a|.")
    print("Given the constraint |a| <= 0.5, the maximum value of 1 + |a| occurs when |a| is maximum, i.e., |a| = 0.5.")
    print("Therefore, the maximum value is 1 + 0.5 = 1.5.")

    print("\nStep 7: Example Polynomial.")
    print("An example polynomial that achieves this maximum is for a = -0.5.")
    a = -0.5
    b = 1.0
    c = 0.5
    print(f"This gives the equation: f(x) = ({a})x^2 + ({b})x + ({c})")
    print(f"For this polynomial, |b| + |c| = |{b}| + |{c}| = {abs(b) + abs(c)}, which is our maximum value.")
    
    max_value = 1.5
    print(f"\nThe final answer for the maximum value of |b| + |c| is {max_value}.")

solve()