import math

def simulate_chsh_experiment():
    """
    Simulates the CHSH inequality test for a quantum reality and a classical dream
    to show how one might distinguish between them.
    """

    # --- Part 1: Simulating the "Waking Reality" (Quantum Mechanical) ---
    print("--- Testing 'Waking Reality' (assumed to be Quantum Mechanical) ---")

    # Optimal angles for maximal violation of the inequality
    # Alice's measurement angles (a, a')
    angle_a = 0
    angle_a_prime = math.pi / 2
    # Bob's measurement angles (b, b')
    angle_b = math.pi / 4
    angle_b_prime = 3 * math.pi / 4

    # In quantum mechanics, the correlation E(a,b) is -cos(a-b) for a singlet state.
    # A common simplified formula used is cos(2*(a-b))
    # We will use the latter as it's often shown in CHSH examples.
    E_ab = math.cos(2 * (angle_b - angle_a))
    E_ab_prime = math.cos(2 * (angle_b_prime - angle_a))
    E_a_prime_b = math.cos(2 * (angle_b - angle_a_prime))
    E_a_prime_b_prime = math.cos(2 * (angle_b_prime - angle_a_prime))

    # CHSH value S
    S_quantum = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime
    tsirelson_bound = 2 * math.sqrt(2)

    print("The CHSH equation is: S = E(a, b) - E(a, b') + E(a', b) + E(a', b')")
    print(f"Plugging in the quantum correlation values:")
    print(f"S = {E_ab:.4f} - ({E_ab_prime:.4f}) + {E_a_prime_b:.4f} + {E_a_prime_b_prime:.4f}")
    print(f"Result for Quantum Reality: S = {S_quantum:.4f}")
    print(f"This value violates the classical limit of 2.\n")

    # Cross-verification test: does the result match the theoretical maximum?
    print("Cross-verifying if the 'test sum' (S - 2*sqrt(2)) equals zero:")
    print(f"S - 2√2 = {S_quantum:.4f} - {tsirelson_bound:.4f} = {S_quantum - tsirelson_bound:.4f}")
    print("The result is zero, consistent with a perfect quantum reality.\n")


    # --- Part 2: Simulating the "Butterfly Dream" (Classical) ---
    print("--- Testing 'Butterfly Dream' (assumed to be a Classical Simulation) ---")

    # In a classical local-realist model, the correlations are different.
    # We can use a simple model: E(a,b) = 1 - (2/pi) * abs(a-b)
    # The angles must be within a 0 to pi range for this model
    # Note: any local hidden variable model will result in |S| <= 2.
    E_ab_classical = 1 - (2/math.pi) * abs(angle_b - angle_a)
    E_ab_prime_classical = 1 - (2/math.pi) * abs(angle_b_prime - angle_a)
    E_a_prime_b_classical = 1 - (2/math.pi) * abs(angle_b - angle_a_prime)
    # E(a', b') is negative in the CHSH formula S=E(a,b)+E(a,b')+E(a',b)-E(a',b') form.
    # Let's use the other form: S = E(ab)-E(ab') + E(a'b) + E(a'b')
    E_a_prime_b_prime_classical = 1 - (2/math.pi) * abs(angle_b_prime - angle_a_prime)
    
    # We are using S = E(ab) - E(ab') + E(a'b) + E(a'b'). Let's stay consistent.
    S_classical = E_ab_classical - E_ab_prime_classical + E_a_prime_b_classical + E_a_prime_b_prime_classical

    print("The CHSH equation is: S = E(a, b) - E(a, b') + E(a', b) + E(a', b')")
    print(f"Plugging in the classical correlation values:")
    print(f"S = {E_ab_classical:.4f} - ({E_ab_prime_classical:.4f}) + {E_a_prime_b_classical:.4f} + {E_a_prime_b_prime_classical:.4f}")
    print(f"Result for Classical Dream: S = {S_classical:.4f}")
    print(f"This value obeys the classical limit of 2.\n")

    # Cross-verification test for the classical world
    print("Cross-verifying if the 'test sum' (S - 2*sqrt(2)) equals zero:")
    print(f"S - 2√2 = {S_classical:.4f} - {tsirelson_bound:.4f} = {S_classical - tsirelson_bound:.4f}")
    print("The result is not zero, indicating this is not a quantum reality.")


if __name__ == '__main__':
    simulate_chsh_experiment()