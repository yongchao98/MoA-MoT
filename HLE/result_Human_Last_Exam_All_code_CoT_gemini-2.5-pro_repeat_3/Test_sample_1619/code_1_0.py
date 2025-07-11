import math
import random

def simulate_chsh_experiment(is_quantum_reality):
    """
    Simulates the result of a CHSH Bell test to determine the nature of a reality.
    - In a classical reality (like a dream), the result S must be <= 2.
    - In a quantum reality (the physical world), the result S can be > 2.

    Args:
        is_quantum_reality (bool): True if simulating quantum reality, False for classical.

    Returns:
        A tuple containing the four correlation values and the final S value.
    """
    if is_quantum_reality:
        # These are ideal correlation values for maximum violation in Quantum Mechanics
        # S can reach up to 2 * sqrt(2) ~= 2.828
        print("...Reality appears to be governed by non-local quantum principles.")
        E1 = 1 / math.sqrt(2)   # E(a, b)
        E2 = -1 / math.sqrt(2)  # E(a, b')
        E3 = 1 / math.sqrt(2)   # E(a', b)
        E4 = 1 / math.sqrt(2)   # E(a', b')
        # S = E(a, b) - E(a, b') + E(a', b) + E(a', b')
        S = E1 - E2 + E3 + E4
        # Add a tiny bit of random "experimental noise"
        S -= random.uniform(0, 0.05)
        return E1, E2, E3, E4, S
    else: # is_classical_reality (dream)
        # These correlation values will respect the classical limit (S <= 2)
        print("...Reality appears to be governed by local, classical principles.")
        E1 = 0.5
        E2 = -0.5
        E3 = 0.5
        E4 = 0.5
        # S = E(a, b) - E(a, b') + E(a', b) + E(a', b')
        S = E1 - E2 + E3 + E4
        # Add a tiny bit of random "experimental noise"
        S -= random.uniform(0, 0.05)
        return E1, E2, E3, E4, S

# --- Main execution ---
# We assume the 'Man' state is the true, quantum reality and the 'Butterfly' is the dream.

# Test 1: The 'Man' Reality
print("--- Testing 'Man' Reality ---")
# Let's assume the Man is in the true, quantum reality.
E1_man, E2_man, E3_man, E4_man, S_man = simulate_chsh_experiment(is_quantum_reality=True)
print("The CHSH equation is: S = E1 - E2 + E3 + E4")
print(f"Final Equation: S = {E1_man:.4f} - ({E2_man:.4f}) + {E3_man:.4f} + {E4_man:.4f} = {S_man:.4f}")
if S_man > 2:
    print(f"Result: S > 2. This violates the Bell inequality, indicating this is the objective, quantum reality.\n")
else:
    print(f"Result: S <= 2. This obeys the Bell inequality, indicating this is a classical simulation (a dream).\n")

# Test 2: The 'Butterfly' Dream
print("--- Testing 'Butterfly' Reality ---")
# Let's assume the Butterfly is in a classical dream state.
E1_bfly, E2_bfly, E3_bfly, E4_bfly, S_bfly = simulate_chsh_experiment(is_quantum_reality=False)
print("The CHSH equation is: S = E1 - E2 + E3 + E4")
print(f"Final Equation: S = {E1_bfly:.4f} - ({E2_bfly:.4f}) + {E3_bfly:.4f} + {E4_bfly:.4f} = {S_bfly:.4f}")
if S_bfly > 2:
    print(f"Result: S > 2. This violates the Bell inequality, indicating this is the objective, quantum reality.")
else:
    print(f"Result: S <= 2. This obeys the Bell inequality, indicating this is a classical simulation (a dream).")
