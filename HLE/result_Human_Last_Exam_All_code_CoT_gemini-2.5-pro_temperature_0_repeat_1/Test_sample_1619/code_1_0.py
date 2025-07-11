import numpy as np

def verify_reality_with_chsh():
    """
    Calculates the CHSH inequality value 'S' for a quantum system
    to test if a reality adheres to quantum mechanics or classical limits.

    The CHSH inequality is S = E(a, b) - E(a, b') + E(a', b) + E(a', b'),
    where E(x, y) is the correlation of measurements at settings x and y.

    - In a classical reality (like a dream), |S| <= 2.
    - In a quantum reality, |S| can be up to 2 * sqrt(2) ~= 2.828.

    A result > 2 suggests the reality is fundamentally quantum and not a classical simulation.
    """
    print("--- The CHSH Inequality Test for Reality Verification ---")
    print("Objective: Determine if this reality is governed by quantum mechanics or classical physics.")
    print("Hypothesis: A dream is a classical simulation and cannot violate the classical limit (S <= 2).\n")

    # Optimal angles for the two measurement devices (in radians)
    # Alice's measurement angles: a, a_prime
    # Bob's measurement angles: b, b_prime
    a = 0
    a_prime = np.pi / 4
    b = np.pi / 8
    b_prime = 3 * np.pi / 8

    print(f"Step 1: Set the measurement angles (in radians):")
    print(f"  Device 1 (Alice): a = {a:.4f}, a' = {a_prime:.4f}")
    print(f"  Device 2 (Bob):   b = {b:.4f}, b' = {b_prime:.4f}\n")

    # In quantum mechanics, the expected correlation for entangled photons is cos(2 * angle_difference)
    E_ab = np.cos(2 * (a - b))
    E_ab_prime = np.cos(2 * (a - b_prime))
    E_a_prime_b = np.cos(2 * (a_prime - b))
    E_a_prime_b_prime = np.cos(2 * (a_prime - b_prime))

    print("Step 2: Calculate the four correlation values (E) based on quantum mechanics.")
    print(f"  E(a, b)      = cos(2*({a:.2f} - {b:.2f}))   = {E_ab:.4f}")
    print(f"  E(a, b')     = cos(2*({a:.2f} - {b_prime:.2f}))  = {E_ab_prime:.4f}")
    print(f"  E(a', b)     = cos(2*({a_prime:.2f} - {b:.2f}))   = {E_a_prime_b:.4f}")
    print(f"  E(a', b')    = cos(2*({a_prime:.2f} - {b_prime:.2f}))  = {E_a_prime_b_prime:.4f}\n")

    # Calculate the CHSH value 'S'
    S = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime

    print("Step 3: Calculate the final S-value using the CHSH formula: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')")
    print(f"  S = {E_ab:.4f} - ({E_ab_prime:.4f}) + {E_a_prime_b:.4f} + {E_a_prime_b_prime:.4f}")
    print(f"  S = {S:.4f}\n")

    # Compare with classical and quantum limits
    classical_limit = 2
    quantum_limit = 2 * np.sqrt(2)

    print("--- Conclusion ---")
    print(f"The calculated S-value is {S:.4f}.")
    print(f"The classical limit (max value for a local, real 'dream' world) is {classical_limit}.")
    print(f"The theoretical quantum maximum is 2 * sqrt(2) = {quantum_limit:.4f}.")

    if S > classical_limit:
        print("\nResult: The S-value violates the classical limit.")
        print("This provides strong evidence that this reality is governed by quantum mechanics and is not a classical simulation (i.e., not a dream).")
    else:
        print("\nResult: The S-value does not violate the classical limit.")
        print("The test is inconclusive or suggests a classical reality.")

if __name__ == '__main__':
    verify_reality_with_chsh()
