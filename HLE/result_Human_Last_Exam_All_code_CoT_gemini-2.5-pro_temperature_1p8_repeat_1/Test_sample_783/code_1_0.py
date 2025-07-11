import math

# Step 1: Define the public parameters
p = 1009  # Prime number
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

print(f"Publicly known information:")
print(f"Prime (p): {p}")
print(f"Generator (g): {g}")
print(f"Alice's public value (A): {A}")
print(f"Bob's public value (B): {B}")
print("-" * 30)

# Step 2: Find Alice's secret number 'a' by brute-force
# We are solving for 'a' in the equation: g^a mod p = A
a = None
print(f"Searching for Alice's secret number 'a'...")
# Iterate from 1 up to p-1
for i in range(1, p):
    # Calculate g^i mod p and check if it equals A
    if pow(g, i, p) == A:
        a = i
        print(f"Found Alice's secret number: a = {a}")
        break

# Step 3: Calculate the shared secret key
if a is not None:
    # The shared secret key S = B^a mod p
    secret_key = pow(B, a, p)
    print("Calculating the shared secret key (S = B^a mod p)...")
    
    # Print the final equation with all the numbers
    print(f"S = {B}^{a} mod {p} = {secret_key}")
    
    print("-" * 30)
    print(f"The final secret key is: {secret_key}")
else:
    print("Could not find Alice's secret number 'a'.")

<<<111>>>