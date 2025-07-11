import random

def run_simulation(N, p0, s, generations):
    """
    Simulates allele frequency change using a simple Wright-Fisher model.
    
    Args:
        N (int): Population size.
        p0 (float): Initial allele frequency.
        s (float): Selection coefficient (s=0 for neutral, s>0 for positive selection).
        generations (int): Number of generations to simulate.
        
    Returns:
        list: A list of allele frequencies over time.
    """
    p = p0
    frequencies = [p]
    
    for _ in range(generations):
        # If allele frequency hits 0 or 1, it's fixed.
        if p == 0 or p == 1:
            frequencies.append(p)
            continue
            
        # Step 1: Calculate frequency after selection
        # The relative fitnesses are 1 for the old allele and 1+s for the new one.
        w_avg = p * (1 + s) + (1 - p) * 1
        p_selected = (p * (1 + s)) / w_avg
        
        # Step 2: Apply genetic drift
        # The number of new alleles in the next generation is a binomial draw.
        # The population has 2N gene copies for a diploid organism.
        num_new_alleles = random.binomial(2 * N, p_selected)
        p = num_new_alleles / (2 * N)
        
        frequencies.append(p)
        
    return frequencies

# Simulation Parameters
N = 500         # Population size
p0 = 0.01       # Initial frequency of the new allele
generations = 100 # Number of generations

# --- Scenario 1: Drift Only (Neutral Allele) ---
s_drift = 0.0
frequencies_drift = run_simulation(N, p0, s_drift, generations)
final_p_drift = frequencies_drift[-1]

# --- Scenario 2: Drift + Positive Selection (Advantageous Allele) ---
s_selection = 0.05 # 5% fitness advantage
frequencies_selection = run_simulation(N, p0, s_selection, generations)
final_p_selection = frequencies_selection[-1]

print("This simulation demonstrates why adaptive evolution challenges models of genetic drift.")
print("-" * 70)
print(f"Parameters: Population Size = {N}, Generations = {generations}, Initial Allele Frequency = {p0}\n")

# Print results for Scenario 1
print(f"Scenario 1: 'Drift Only'")
print(f"The selection coefficient is: s = {s_drift}")
print(f"With only random drift acting, the final allele frequency after {generations} generations is: {final_p_drift:.4f}")
print("Under pure drift, a new allele is very likely to be lost by random chance.\n")

# Print results for Scenario 2
print(f"Scenario 2: 'Drift + Positive Selection'")
print(f"The selection coefficient is: s = {s_selection}")
print(f"When positive selection is added, the final allele frequency after {generations} generations is: {final_p_selection:.4f}")
print("Adaptive evolution provides a strong directional force, 'outweighing' drift and driving the allele to a higher frequency.")
print("-" * 70)

<<<C>>>