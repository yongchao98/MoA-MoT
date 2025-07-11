# Demonstrating the truth of statement G.
# G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions

# Let's choose some positive integers for n and m.
n = 17
m = 5

# Perform the floor division and modulo operations.
quotient = n // m
remainder = n % m

# Reconstruct n using the division identity.
reconstructed_n = quotient * m + remainder

# Print out the step-by-step calculation.
print(f"Verifying the identity (n//m)*m + n%m == n for n={n} and m={m}:")
print(f"1. Start with the expression: ({n}//{m})*{m} + {n}%{m}")
print(f"2. Calculate the quotient and remainder: {quotient} and {remainder}")
print(f"3. Substitute them into the expression: ({quotient})*{m} + {remainder}")
print(f"4. Perform the multiplication: {quotient * m} + {remainder}")
print(f"5. Perform the addition: {reconstructed_n}")
print(f"6. Compare with the original n: Is {reconstructed_n} == {n}?")

# Final verification
is_true = reconstructed_n == n
print(f"The statement holds: {is_true}")