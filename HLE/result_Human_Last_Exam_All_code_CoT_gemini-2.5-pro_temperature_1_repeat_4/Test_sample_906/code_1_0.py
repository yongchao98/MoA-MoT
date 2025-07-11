import math

# Define rho symbolically for the formula
rho_symbol = "\u03C1" # Unicode for Greek small letter rho

# The derivation steps lead to the final equation for the steady-state probability pi_0.
# 1. The balance equations give pi_k = pi_0 * (rho^k / k!).
# 2. The normalization condition sum(pi_k) = 1 gives pi_0 * sum(rho^k / k!) = 1.
# 3. The sum is the Taylor series for e^rho.
# 4. Therefore, pi_0 * e^rho = 1.

# Final Equation
pi_0 = f"e^(-{rho_symbol})"

# Print the final result
print(f"The steady-state probability pi_0 in terms of {rho_symbol} = \u03BB/\u03BC is:")
print(f"\u03C0_0 = {pi_0}")