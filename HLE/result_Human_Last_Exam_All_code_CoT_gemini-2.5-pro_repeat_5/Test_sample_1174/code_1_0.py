# Values of the exponents derived from physical principles.

# n1, n2, n3 are from the zero-frequency limit: S_B ~ T * sigma * t / z^4
# The problem statement implies S_B ~ sigma^n1 * T^n2 * z^n3
n1 = 1
n2 = 1
n3 = -4

# n4, n5, n6 are from the frequency scaling in different regimes.
# Regime 1: omega << 1/(sigma*z*t) -> Unscreened, white noise spectrum
n4 = 0
# Regime 2: 1/(sigma*z*t) << omega << 1/(sigma*t^2) -> Screened thin film
n5 = -2
# Regime 3: omega >> 1/(sigma*t^2) -> Thick film (skin depth limited)
n6 = -1

# The list of exponents
n = [n1, n2, n3, n4, n5, n6]

# Calculate the sum sum_{k=1 to 6} k * n_k
total_sum = 0
equation_parts = []
for k in range(1, 7):
    term = k * n[k-1]
    total_sum += term
    
    # Format the string for each term in the sum
    if n[k-1] < 0:
        equation_parts.append(f"{k}*({n[k-1]})")
    else:
        equation_parts.append(f"{k}*{n[k-1]}")

# Print the full equation
equation_str = " + ".join(equation_parts)
print(f"The calculation is: {equation_str}")

# Print the final result
print(f"The value of the sum is: {total_sum}")
