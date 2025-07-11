import math

def falling_factorial(x, k):
    """Computes the falling factorial (x)_k = x(x-1)...(x-k+1)."""
    if k == 0:
        return 1.0
    res = 1.0
    for i in range(k):
        res *= (x - i)
    return res

def binomial_coefficient(n, k):
    """Computes the binomial coefficient n choose k."""
    if k < 0 or k > n:
        return 0
    # Using math.comb for precision and efficiency
    return math.comb(n, k)

def calculate_Sn(n):
    """Calculates the sum S_n."""
    total_sum = 0
    x = n + 0.5
    for m in range(n + 1):
        term = ((-1)**m) * binomial_coefficient(n, m) * falling_factorial(x, m)
        total_sum += term
    return total_sum

# Calculate and print S_n for n from 0 to 7
for n in range(8):
    sn_value = calculate_Sn(n)
    print(f"n = {n}, S_n = {sn_value}")

# Analyze the growth by comparing |S_n| with simple functions
print("\n--- Growth Analysis ---")
for n in range(1, 8):
    sn_value = calculate_Sn(n)
    abs_sn = abs(sn_value)
    
    # Check against f(n) = 1 (constant)
    is_bounded_by_C = "Yes" if abs_sn < 4 else f"No, value is {abs_sn:.2f}"

    # Check against f(n) = n
    if n > 0:
        ratio_n = abs_sn / n
        print(f"n = {n}, |S_n|/n = {ratio_n:.2f}")

# The problem asks for the function f, not the constant C.
# From the analysis, we can see that the sequence S_n is not bounded by a constant.
# Let's observe the ratio of consecutive terms |S_n / S_{n-1}| to check for exponential growth.
print("\n--- Ratio Analysis ---")
s_prev = calculate_Sn(2) # Start from n=3 as S_0, S_1, S_2 are small
for n in range(3, 8):
    s_curr = calculate_Sn(n)
    if abs(s_prev) > 1e-9: # Avoid division by zero
        ratio = abs(s_curr / s_prev)
        print(f"n = {n}, |S_{n}/S_{n-1}| = {ratio:.2f}")
    s_prev = s_curr