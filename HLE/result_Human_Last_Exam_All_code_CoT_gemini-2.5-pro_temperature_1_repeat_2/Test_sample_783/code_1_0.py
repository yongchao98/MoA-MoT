import math

# Publicly known numbers
p = 1009  # Prime number
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

a_secret = None

# Step 1: Find Alice's secret number 'a' by solving g^a mod p = A
# We will iterate through possible values for 'a' from 1 to p-1.
print(f"Searching for Alice's secret number 'a' such that {g}^a mod {p} = {A}...")

for i in range(1, p):
    # Calculate g^i mod p and check if it matches Alice's public number A
    if pow(g, i, p) == A:
        a_secret = i
        print(f"Found Alice's secret number: a = {a_secret}")
        break

# Step 2: If 'a' was found, calculate the shared secret key S = B^a mod p
if a_secret is not None:
    # Calculate the shared secret using Bob's public number B and Alice's secret 'a'
    shared_secret = pow(B, a_secret, p)
    
    print("\nCalculating the shared secret key...")
    print(f"The final equation is: {B} ^ {a_secret} mod {p} = {shared_secret}")
    print(f"\nThe secret key is: {shared_secret}")
else:
    print("Could not determine Alice's secret number 'a'.")

<<<182>>>