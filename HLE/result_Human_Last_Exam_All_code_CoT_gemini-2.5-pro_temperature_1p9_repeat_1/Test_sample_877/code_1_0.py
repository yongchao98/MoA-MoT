import math

def solve():
    """
    This function determines the function h(x) based on the analysis of the ODE system.
    The condition for a(t) not escaping to -infinity is that a'(0) must be non-negative.
    a'(0) = -1/2*a(0)^2 - 2*b(0)^2 - b(0) + 1 >= 0
    a(0)^2 <= 2 - 2*b(0) - 4*b(0)^2
    The problem asks for a condition of the form -sqrt(h(b(0))) < a(0) < 0,
    which implies a(0)^2 < h(b(0)).
    Thus, we identify h(x) = 2 - 2*x - 4*x^2.
    """
    
    # Coefficients of the quadratic function h(x) = c_0 + c_1*x + c_2*x^2
    c_0 = 2
    c_1 = -2
    c_2 = -4
    
    # We output the function h(x) with its coefficients.
    print("The function h(x) is determined from the condition a'(0) >= 0.")
    print("This leads to the inequality a(0)^2 <= 2 - 2*b(0) - 4*b(0)^2.")
    print("By comparing with the required form -sqrt(h(b(0))) < a(0) < 0, we identify h(x).")
    print("\nThe function h(x) is:")
    print(f"h(x) = {c_0} + ({c_1})*x + ({c_2})*x^2")
    
    final_h_expression = f"{c_0} - {abs(c_1)}*x - {abs(c_2)}*x**2"
    # Provide the answer in the specified format
    print(f"\nFinal answer in formula form: <<<h(x) = {final_h_expression}>>>")

solve()

# The final answer content required is the mathematical expression for h(x)
# Extracted from the code's logic.
# h(x) = 2 - 2x - 4x^2