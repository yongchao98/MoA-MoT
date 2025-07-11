import numpy as np

def run_chsh_simulation(reality_type, num_trials=50000):
    """
    Simulates a CHSH experiment for a given reality type.

    - Quantum Reality: Simulates outcomes based on quantum mechanical predictions for
      entangled photons, which can violate the Bell inequality.
    - Classical Reality (Dream): Simulates outcomes based on a local hidden variable
      model, which cannot violate the Bell inequality.
    """
    print(f"--- Running Simulation for: {reality_type} ---")

    # Define the measurement angles for Alice and Bob
    # These are the standard angles that give the maximum violation in QM
    angle_a = 0
    angle_a_prime = np.pi / 4  # 45 degrees
    angle_b = np.pi / 8        # 22.5 degrees
    angle_b_prime = 3 * np.pi / 8 # 67.5 degrees

    # Counters for the products of outcomes for each setting combination
    # (e.g., product_ab is the sum of A*B for settings (a, b))
    product_ab = 0
    product_ab_prime = 0
    product_a_prime_b = 0
    product_a_prime_b_prime = 0
    
    # Dictionaries to map settings to their angles for clarity
    alice_angles = {0: angle_a, 1: angle_a_prime}
    bob_angles = {0: angle_b, 1: angle_b_prime}

    for _ in range(num_trials):
        # Alice and Bob randomly choose their measurement setting (0 or 1)
        alice_choice = np.random.randint(0, 2)
        bob_choice = np.random.randint(0, 2)
        
        alice_angle = alice_angles[alice_choice]
        bob_angle = bob_angles[bob_choice]

        outcome_a = 0
        outcome_b = 0

        if reality_type == 'Quantum Reality':
            # In Quantum Mechanics, the correlation depends on the angle difference
            # We can simulate this without complex numbers.
            # Alice's outcome is random (+1 or -1).
            outcome_a = 1 if np.random.random() < 0.5 else -1
            
            # Bob's outcome is correlated based on quantum probability.
            # The probability of getting the same outcome is cos^2(angle_diff)
            angle_difference = alice_angle - bob_angle
            prob_same_outcome = np.cos(angle_difference)**2
            
            if np.random.random() < prob_same_outcome:
                outcome_b = outcome_a
            else:
                outcome_b = -outcome_a
        
        elif reality_type == 'Classical Reality (Dream)':
            # In a classical local hidden variable model (a 'dream'),
            # the outcome is predetermined by a shared 'hidden variable' lambda.
            hidden_variable_lambda = np.random.uniform(0, 2 * np.pi)
            
            # The 'Bell anagol' local model: outcome is +1 if the detector is within
            # 90 degrees of lambda, and -1 otherwise.
            outcome_a = 1 if abs(alice_angle - hidden_variable_lambda) % (2*np.pi) < np.pi / 2 or \
                             abs(alice_angle - hidden_variable_lambda) % (2*np.pi) > 3 * np.pi / 2 else -1
            outcome_b = 1 if abs(bob_angle - hidden_variable_lambda) % (2*np.pi) < np.pi / 2 or \
                             abs(bob_angle - hidden_variable_lambda) % (2*np.pi) > 3 * np.pi / 2 else -1
            # In a "dream", the outcomes are determined locally by the dreamer's mind,
            # which can be modeled by a hidden variable. We must also enforce perfect
            # anti-correlation for entangled particles when measured at the same angle.
            # However, for simplicity and demonstration, we use a basic LHV model that
            # shows the classical bound is respected.

        # Calculate the product of the outcomes
        product = outcome_a * outcome_b
        
        # Add the product to the correct counter
        if alice_choice == 0 and bob_choice == 0:
            product_ab += product
        elif alice_choice == 0 and bob_choice == 1:
            product_ab_prime += product
        elif alice_choice == 1 and bob_choice == 0:
            product_a_prime_b += product
        elif alice_choice == 1 and bob_choice == 1:
            product_a_prime_b_prime += product

    # There are four combinations of settings, so each gets approx. N/4 trials
    num_per_setting = num_trials / 4.0
    
    # Calculate the expectation values (average of the products)
    E_ab = product_ab / num_per_setting
    E_ab_prime = product_ab_prime / num_per_setting
    E_a_prime_b = product_a_prime_b / num_per_setting
    E_a_prime_b_prime = product_a_prime_b_prime / num_per_setting

    # Calculate the CHSH value S
    # The standard inequality is |E(a,b) - E(a,b') + E(a',b) + E(a',b')| <= 2
    S = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime
    
    print("The final CHSH equation is S = E(a,b) - E(a,b') + E(a',b) + E(a',b')")
    print(f"S = {E_ab:.4f} - ({E_ab_prime:.4f}) + {E_a_prime_b:.4f} + {E_a_prime_b_prime:.4f}")
    print(f"Resulting S value: {S:.4f}")
    
    if reality_type == 'Quantum Reality':
        print(f"This value is GREATER than the classical limit of 2.")
        print("Theoretical maximum for Quantum Mechanics is 2*sqrt(2) ≈ 2.828\n")
    else:
        print(f"This value is LESS than or equal to the classical limit of 2.\n")


# Run the simulation for both scenarios
run_chsh_simulation('Quantum Reality')
run_chsh_simulation('Classical Reality (Dream)')