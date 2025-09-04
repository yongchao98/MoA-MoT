import random
import numpy as np

def solve_magnetic_energy_comparison():
    """
    Uses Monte Carlo sampling and deterministic verification to compare
    paramagnetic coupling energy with a hydrogen transition energy.
    """
    # --- Setup: Define constants and parameters ---
    # Physical constants
    mu_B_eV = 5.7883818012e-5  # Bohr magneton in eV/T
    hc_eV_nm = 1239.84198      # Planck's constant * speed of light in eV*nm
    
    # Given parameters
    B = 1.0  # Magnetic field in Tesla
    lambda_um = 0.4861 # Wavelength in micrometers
    lambda_nm = lambda_um * 1000 # Convert wavelength to nanometers

    # --- (a) Sampling / Monte Carlo Exploration ---
    num_trials = 100
    # The problem states "small values of m". We'll sample from a representative range.
    # The transition lambda=486.1nm is from n=4 to n=2. For n=4, l can be up to 3,
    # so m can be up to 3. We'll sample m from {1, 2, 3}.
    possible_m = range(1, 4) 
    
    print("--- (a) Sampling / Monte Carlo Exploration ---")
    print(f"Running {num_trials} trials, sampling 'm' from {list(possible_m)}...")
    
    # Store comparison results
    # Options are: A) =, B) >>, C) <<, D) >
    candidate_votes = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    
    # Calculate the constant transition energy once
    delta_E = hc_eV_nm / lambda_nm

    for _ in range(num_trials):
        m = random.choice(possible_m)
        
        # Calculate paramagnetic coupling energy for the sampled m
        H_para = m * mu_B_eV * B
        
        # Compare magnitudes
        ratio = H_para / delta_E
        
        if ratio < 1e-2: # Threshold for "much less than" (<<)
            candidate_votes['C'] += 1
        elif ratio > 1e2: # Threshold for "much greater than" (>>)
            candidate_votes['B'] += 1
        elif np.isclose(ratio, 1.0): # Threshold for "equal to" (=)
            candidate_votes['A'] += 1
        elif ratio > 1.0: # Greater than, but not "much greater" (>)
            candidate_votes['D'] += 1
        # Note: The case for '<' is covered by '<<' as the ratio is extremely small.

    print(f"Sampling results (votes for each option): {candidate_votes}")

    # --- (b) Narrowing Candidates ---
    print("\n--- (b) Narrowing Candidates ---")
    if not any(candidate_votes.values()):
        most_likely_option = None
        print("No clear candidate emerged from sampling.")
    else:
        most_likely_option = max(candidate_votes, key=candidate_votes.get)
        print(f"The most likely candidate from sampling is Option {most_likely_option}.")

    # --- (c) Deterministic Verification ---
    print("\n--- (c) Deterministic Verification ---")
    
    # Use a representative small, non-zero value for m, m=1
    m_det = 1
    
    # Calculate paramagnetic coupling energy <H>
    H_para_det = m_det * mu_B_eV * B
    
    # The transition energy Delta E is the same as before
    delta_E_det = delta_E
    
    # Calculate the precise ratio
    final_ratio = H_para_det / delta_E_det
    
    print(f"Using representative value m = {m_det}:")
    print(f"Paramagnetic coupling energy <H> = {H_para_det:.4e} eV")
    print(f"Transition energy Delta E = {delta_E_det:.4f} eV")
    print(f"Ratio <H> / Delta E = {final_ratio:.4e}")
    
    # Verify the relationship against the options
    verified_answer = None
    if final_ratio < 1e-2:
        verified_relationship = "<< (much less than)"
        verified_answer = "C"
    elif final_ratio > 1e2:
        verified_relationship = ">> (much greater than)"
        verified_answer = "B"
    elif np.isclose(final_ratio, 1.0):
        verified_relationship = "= (equal to)"
        verified_answer = "A"
    elif final_ratio > 1:
        verified_relationship = "> (greater than)"
        verified_answer = "D"
    else:
        verified_relationship = "undetermined"

    # --- (d) Print Final Verified Answer ---
    print("\n--- (d) Final Verified Answer ---")
    if verified_answer and verified_answer == most_likely_option:
        print(f"The deterministic check confirms the sampling result.")
        print(f"The paramagnetic energy ({H_para_det:.2e} eV) is several orders of magnitude smaller than the transition energy ({delta_E_det:.2f} eV).")
        print(f"The relationship is <H> {verified_relationship} Delta E.")
        print(f"The correct option is {verified_answer}.")
    elif not verified_answer:
        print("Could not deterministically verify any of the given options.")
    else:
        print(f"Verification ({verified_answer}) contradicts sampling ({most_likely_option}). Please check the logic.")

# Execute the full analysis
solve_magnetic_energy_comparison()