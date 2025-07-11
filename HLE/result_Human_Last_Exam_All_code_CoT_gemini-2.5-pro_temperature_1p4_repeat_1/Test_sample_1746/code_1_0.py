# This script verifies statement (G) using an example.
# Statement G: For positive integers n, m: (n//m)*m + n%m == n holds true in all versions.

# 1. Define positive integers for the dividend (n) and divisor (m).
n = 27
m = 4

# 2. Calculate the quotient (q) and remainder (r) as Python does.
q = n // m
r = n % m

# 3. Reconstruct the dividend using the calculated quotient and remainder.
reconstructed_n = q * m + r

# 4. Print the verification process step-by-step, showing all numbers involved.
print(f"Verifying the identity: (n // m) * m + n % m == n")
print(f"With n = {n} and m = {m}:")
print(f"({n} // {m}) * {m} + ({n} % {m}) == {n}")
print("")
print(f"First, we evaluate the quotient and remainder:")
print(f"  Quotient ({n} // {m}) is {q}")
print(f"  Remainder ({n} % {m}) is {r}")
print("")
print(f"Next, we substitute these values into the equation:")
print(f"({q}) * {m} + ({r}) == {n}")
print("")
print(f"Finally, we evaluate the left side of the equation:")
print(f"{q * m} + {r} == {n}")
print(f"{reconstructed_n} == {n}")
print("")

# 5. Conclude with the result of the boolean comparison.
print(f"Is the identity true for these numbers? {reconstructed_n == n}")
