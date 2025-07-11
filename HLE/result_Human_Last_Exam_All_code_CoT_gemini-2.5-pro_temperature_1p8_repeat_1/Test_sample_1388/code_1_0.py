import math

# The derivation leads to the conclusion that H(t) is an exponential function
# of h(t), which is the integrated L-infinity norm of the spatial derivative of u.
# h(t) = integral from 0 to t of ||∂_x u(·, s)||_L∞ ds

# The final form of H(t) is determined as:
H_t_expression = "exp(h(t))"
H_t_full_expression = "exp(Integral from 0 to t of ||∂_x u(·, s)||_L∞ ds)"

print("Based on the L2 energy estimate, the explicit form of H(t) is:")
print(f"H(t) = {H_t_expression}")
print("Which can be written in full as:")
print(f"H(t) = {H_t_full_expression}")

# The provided inequality is:
# ||u(·, t)||_L2 ≤ ||u_0||_L2 * H(t)
# We found:
# ||u(·, t)||_L2 ≤ ||u_0||_L2 * exp(h(t))
# Therefore, H(t) is exp(h(t)).

final_answer = "exp(integral_from_0_to_t ||u_x(s)||_L_infinity ds)"
# There are no numerical constants other than the implicit '1' and '0' in the expression.
# So we present the expression itself as the output.
print("\nThe final result for H(t) is:")
# Printing each component of the expression as requested in the instructions
# "Remember in the final code you still need to output each number in the final equation!"
# The expression is symbolic, we will represent it in text.
print("H(t) = exp(integral from 0 to t of ||∂_x u(s)||_L_infinity ds)")
