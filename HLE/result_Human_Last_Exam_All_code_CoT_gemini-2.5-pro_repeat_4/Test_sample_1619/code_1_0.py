import numpy as np

def simulate_chsh_experiment(num_trials_per_setting, reality_type):
    """
    Simulates a CHSH Bell test experiment for a given reality type.

    In this experiment, a source sends entangled particles to two observers,
    Alice and Bob. They choose measurement settings and compare results.

    Alice's settings: a=0°, a'=90°
    Bob's settings:   b=45°, b'=135°

    These angles are chosen to maximally violate the inequality in the quantum case.
    """
    # Define the measurement angles in radians
    alice_settings = [0, np.pi / 2]  # Settings a and a'
    bob_settings = [np.pi / 4, 3 * np.pi / 4]  # Settings b and b'

    # This dictionary will store the calculated correlation for each pair of settings.
    correlations = {}

    # We test all four combinations of settings
    for a_setting in alice_settings:
        for b_setting in bob_settings:
            product_sum = 0
            for _ in range(num_trials_per_setting):
                # In a real experiment, particles are measured one by one.
                # Here, we simulate the outcome for each trial.
                alice_outcome = 0
                bob_outcome = 0

                if reality_type == 'quantum':
                    # In Quantum Mechanics, the correlation is a function of the
                    # difference between the measurement angles.
                    # Correlation E = -cos(angle_difference)
                    angle_diff = a_setting - b_setting
                    quantum_correlation = -np.cos(angle_diff)
                    
                    # We simulate outcomes based on this correlation.
                    # Probability that outcomes have a product of +1 is (1+E)/2
                    if np.random.random() < (1 + quantum_correlation) / 2:
                        product_sum += 1
                    else:
                        product_sum -= 1

                elif reality_type == 'classical':
                    # In a Classical (Local Hidden Variable) model, the particle's
                    # properties are predetermined by a "hidden variable", here an angle.
                    hidden_variable_angle = np.random.uniform(0, 2 * np.pi)
                    
                    # Define a simple deterministic outcome based on the hidden variable.
                    # Outcome is +1 if the setting is within 90° of the hidden angle.
                    def get_classical_outcome(setting, hidden_angle):
                        diff = (setting - hidden_angle + np.pi) % (2 * np.pi) - np.pi
                        return 1 if abs(diff) <= np.pi / 2 else -1

                    alice_outcome = get_classical_outcome(a_setting, hidden_variable_angle)
                    # For perfect anti-correlation, Bob's result for a given setting
                    # must be the opposite of Alice's.
                    bob_outcome = -get_classical_outcome(b_setting, hidden_variable_angle)
                    
                    product_sum += alice_outcome * bob_outcome
            
            # The correlation is the average of the product of outcomes.
            correlations[(a_setting, b_setting)] = product_sum / num_trials_per_setting

    # Extract the four correlation values needed for the CHSH 'S' value
    E_ab = correlations[(alice_settings[0], bob_settings[0])]
    E_ab_prime = correlations[(alice_settings[0], bob_settings[1])]
    E_a_prime_b = correlations[(alice_settings[1], bob_settings[0])]
    E_a_prime_b_prime = correlations[(alice_settings[1], bob_settings[1])]

    # Calculate the CHSH 'S' value using the formula: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
    s_value = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime
    
    return s_value, E_ab, E_ab_prime, E_a_prime_b, E_a_prime_b_prime

# --- Main ---
# Set the number of simulated particle pairs for the experiment
num_trials = 200000

# Simulate the test in the "Quantum" reality
S_quantum, q1, q2, q3, q4 = simulate_chsh_experiment(num_trials, 'quantum')

# Simulate the test in the "Classical" (Dream) reality
S_classical, c1, c2, c3, c4 = simulate_chsh_experiment(num_trials, 'classical')

# --- Print The Results ---
print("--- Test in the 'Man' (Quantum) Reality ---")
print(f"Based on {num_trials} trials per setting, the correlations are:")
print(f"E(a, b)   = {q1:.4f}\nE(a, b')  = {q2:.4f}\nE(a', b)  = {q3:.4f}\nE(a', b') = {q4:.4f}")
print("\nThe CHSH equation is: S = E(a, b) - E(a, b') + E(a', b) + E(a', b')")
print(f"Plugging in the numbers: S = ({q1:.4f}) - ({q2:.4f}) + ({q3:.4f}) + ({q4:.4f})")
print(f"Final S-value: {S_quantum:.4f}")
print(f"Result: The absolute value |S| is {abs(S_quantum):.4f}, which violates the classical limit of 2.\n")

print("="*50 + "\n")

print("--- Test in the 'Butterfly' (Classical Dream) Reality ---")
print(f"Based on {num_trials} trials per setting, the correlations are:")
print(f"E(a, b)   = {c1:.4f}\nE(a, b')  = {c2:.4f}\nE(a', b)  = {c3:.4f}\nE(a', b') = {c4:.4f}")
print("\nThe CHSH equation is: S = E(a, b) - E(a, b') + E(a', b) + E(a', b')")
print(f"Plugging in the numbers: S = ({c1:.4f}) - ({c2:.4f}) + ({c3:.4f}) + ({c4:.4f})")
print(f"Final S-value: {S_classical:.4f}")
print(f"Result: The absolute value |S| is {abs(S_classical):.4f}, which obeys the classical limit of 2.\n")

print("--- Conclusion ---")
print("The test provides a clear distinction. The reality that violates the classical limit is the true quantum reality, while the other is a simulated dream.")
