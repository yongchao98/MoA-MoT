import sympy as sp

def solve():
    """
    This function derives and prints the function h(b) that determines the stability condition.
    The derivation steps are outlined in the text explanation.
    The final result is found to be h(b) = 2 - 2*b - 4*b**2.
    """
    
    # Define the symbol for the variable b
    b = sp.Symbol('b')
    
    # The derived function h(b)
    h_b = 2 - 2*b - 4*b**2
    
    # We want to present the answer as an equation for h(b)
    # The coefficients of the polynomial h(b) = c2*b^2 + c1*b + c0
    # can be extracted if needed, but the prompt asks to output the equation.
    
    coeffs = h_b.as_poly(b).all_coeffs()
    c2, c1, c0 = coeffs[0], coeffs[1], coeffs[2]
    
    print("The function h(b) is found by analyzing the system's Lyapunov function.")
    print("The condition for a(t) to approach 0 is -sqrt(h(b(0))) < a(0) < 0, where h(b) is:")
    print(f"h(b) = {c0} + ({c1})*b + ({c2})*b^2")
    
    # As requested, output each number in the final equation.
    print("\nIn full form:")
    print(f"h(b) = {c0} - {abs(c1)}*b - {abs(c2)}*b^2")

solve()

# The final expression for h(b) is provided below, surrounded by "<<<" and ">>>".
# We represent the function h(x) in terms of a variable x.
# h(x) = 2 - 2*x - 4*x^2
final_expression = "2 - 2*x - 4*x**2"
print(f"\n<<<h(x) = {final_expression}>>>")
