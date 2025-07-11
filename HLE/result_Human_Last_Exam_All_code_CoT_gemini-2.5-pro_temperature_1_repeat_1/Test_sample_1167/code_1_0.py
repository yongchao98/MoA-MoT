import math

# Step 1: Define the problem parameters
# We are given a threshold T = N^(3/8)
# The set is X = {x: max_t |S(x,t)| > T}
# We want to find alpha where |X| <= C * N^alpha

# Step 2: Use Chebyshev's inequality
# |X| * T^p <= Integral(M(x)^p dx) for any p > 0.
# Let's choose p=4.
# |X| * (N^(3/8))^4 <= Integral(M(x)^4 dx)
# |X| * N^(3/2) <= Integral(M(x)^4 dx)

# Step 3: Bound the integral of M(x)^4
# We showed in the plan that Integral(M(x)^4 dx) is bounded by O(N^2).
# Let's state this as an inequality: Integral(M(x)^4 dx) <= K * N^2 for some constant K.

# Step 4: Combine the inequalities to find the bound on |X|
# |X| * N^(3/2) <= K * N^2
# |X| <= K * N^2 / N^(3/2)
# |X| <= K * N^(2 - 3/2)
# |X| <= K * N^(1/2)

# Step 5: Determine alpha
# The upper bound for |X| is of the form N^alpha.
# From our derivation, alpha = 1/2.

# The problem asks for the real number alpha.
alpha_numerator = 1
alpha_denominator = 2
alpha = alpha_numerator / alpha_denominator

# Let's represent the logic in the code
N_power_in_threshold = 3/8
p_for_chebyshev = 4

# The exponent of N in the Chebyshev inequality term is:
N_exponent_chebyshev = N_power_in_threshold * p_for_chebyshev
# 3/8 * 4 = 12/8 = 3/2

# The exponent of N in the moment integral bound is:
# Integral(M(x)^4 dx) = O(N^2), so the exponent is 2.
N_exponent_integral = 2

# The final exponent alpha is calculated as:
# |X| * N^(N_exponent_chebyshev) <= K * N^(N_exponent_integral)
# |X| <= K * N^(N_exponent_integral - N_exponent_chebyshev)
alpha = N_exponent_integral - N_exponent_chebyshev

# Final calculation
final_numerator = 2 * 2 - 3
final_denominator = 2
alpha_val = final_numerator / final_denominator

# We are asked to output the logic and the final number.
# The calculation shows:
# alpha = 2 - (3/8)*4
# alpha = 2 - 12/8
# alpha = 2 - 3/2
# alpha = 4/2 - 3/2
# alpha = 1/2

print("The best upper bound for |X| is of the form N^alpha.")
print("We use Chebyshev's inequality with the 4th moment:")
print("|X| * (N^(3/8))^4 <= integral( (max_t |S(x,t)|)^4 dx )")
print("This simplifies to: |X| * N^(3/2) <= integral( M(x)^4 dx )")
print("Using Nikolskii's inequality and properties of the sum, we can show:")
print("integral( M(x)^4 dx ) <= C * N^2")
print("Combining these gives: |X| * N^(3/2) <= C * N^2")
print("So, |X| <= C * N^(2 - 3/2) = C * N^(1/2)")
print("This implies that the exponent alpha is 1/2.")
print("The value of alpha is determined by the equation: alpha = 2 - 4 * (3/8)")
print(f"alpha = {N_exponent_integral} - {p_for_chebyshev} * (3/8) = {alpha}")
