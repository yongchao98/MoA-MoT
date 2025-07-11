import math

def final_answer():
    """
    This function prints the derived function h(x) and its coefficients.
    The condition for a(t) to converge to 0 is:
    -sqrt(h(b(0))) < a(0) < 0
    """

    # Coefficients and constants from the derivation of h(x)
    # h(x) = c1*x**2 + c2*x + c3*ln(x) + num1 + num2
    c1 = 2.0
    c2 = 2.0
    c3 = -2.0
    num1 = -1.5  # This is -3/2
    num2 = -2.0 * math.log(2.0)
    
    # The problem requires outputting each number in the final equation.
    print("The function h(x) that defines the convergence region is:")
    print(f"h(x) = ({c1})*x^2 + ({c2})*x + ({c3})*ln(x) + ({num1}) + ({num2})")
    print("\nFor clarity, the symbolic form is:")
    print("h(x) = 2*x^2 + 2*x - 2*ln(x) - 3/2 - 2*ln(2)")

final_answer()