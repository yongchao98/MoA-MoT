def solve_diffie_hellman():
    """
    Solves for the shared secret key in a Diffie-Hellman exchange
    by brute-forcing one of the private keys.
    """
    # Given parameters
    p = 1009         # Prime number
    g = 11           # Base
    alice_public = 297 # Alice's public number (g^a mod p)
    bob_public = 944   # Bob's public number (g^b mod p)

    a_secret = None

    # Brute-force to find Alice's secret number 'a'
    # We are looking for 'a' such that: g^a mod p = alice_public
    for i in range(1, p):
        if pow(g, i, p) == alice_public:
            a_secret = i
            break

    if a_secret is not None:
        # Once we have 'a', we can calculate the shared secret key
        # S = (bob_public)^a mod p
        shared_secret = pow(bob_public, a_secret, p)
        
        # Print the final equation with all the numbers
        print("To find the secret key, we first find Alice's secret number 'a'.")
        print(f"By checking powers of g ({g}), we find that a = {a_secret} because {g}^{a_secret} mod {p} = {alice_public}.")
        print("\nThen, we can calculate the shared secret key using Bob's public number:")
        print(f"Equation: {bob_public} ^ {a_secret} mod {p} = {shared_secret}")
        print(f"The shared secret key is: {shared_secret}")
    else:
        print("Could not find Alice's secret number 'a' within the given range.")

if __name__ == '__main__':
    solve_diffie_hellman()
<<<277>>>