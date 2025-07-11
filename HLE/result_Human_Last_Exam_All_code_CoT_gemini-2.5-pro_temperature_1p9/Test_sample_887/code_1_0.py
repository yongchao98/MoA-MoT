import sys
# This script simulates the effect of bond-selective laser excitation
# on the chemical reactivity of CHD3 with atomic fluorine.

def simulate_reaction():
    """
    Models the reaction F + CHD3 -> HF + CD3 or DF + CHD2
    both with and without laser excitation of the C-H bond.
    """
    # --- Parameters ---
    num_ch_bonds = 1
    num_cd_bonds = 3

    # Baseline reactivity. The C-D bond is slightly stronger (and thus less reactive)
    # than the C-H bond due to a lower zero-point energy. We'll use an illustrative value.
    base_reactivity_ch = 1.0
    base_reactivity_cd = 0.8

    # The laser deposits energy into the C-H bond, dramatically increasing its reactivity.
    # Experimental results show an enhancement factor of about 30.
    excitation_factor = 30.0

    print("--- Simulating the reaction of F + CHD3 ---\n")

    # --- Scenario 1: No Laser ---
    print("SCENARIO 1: No Laser Excitation")
    
    # Calculate the total reactivity for H and D atom abstraction channels
    h_abstraction_score_no_laser = num_ch_bonds * base_reactivity_ch
    d_abstraction_score_no_laser = num_cd_bonds * base_reactivity_cd
    total_reactivity_no_laser = h_abstraction_score_no_laser + d_abstraction_score_no_laser
    
    print(f"Total relative H-removal score: {num_ch_bonds} (bond) * {base_reactivity_ch:.1f} (reactivity) = {h_abstraction_score_no_laser:.1f}")
    print(f"Total relative D-removal score: {num_cd_bonds} (bonds) * {base_reactivity_cd:.1f} (reactivity) = {d_abstraction_score_no_laser:.1f}")
    
    # Calculate the probability of removing an H atom
    prob_h_removal_no_laser = h_abstraction_score_no_laser / total_reactivity_no_laser
    print(f"Result: The likelihood of H atom removal is {prob_h_removal_no_laser:.1%}.\n")

    # --- Scenario 2: With Laser Exciting the C-H Bond ---
    print("SCENARIO 2: With Laser Excitation of the C-H Bond")
    
    # Apply the excitation factor to the C-H bond's reactivity
    reactivity_ch_excited = base_reactivity_ch * excitation_factor
    
    h_abstraction_score_laser = num_ch_bonds * reactivity_ch_excited
    d_abstraction_score_laser = num_cd_bonds * base_reactivity_cd # C-D reactivity is unaffected
    total_reactivity_laser = h_abstraction_score_laser + d_abstraction_score_laser

    print(f"Total relative H-removal score: {num_ch_bonds} (bond) * {reactivity_ch_excited:.1f} (reactivity) = {h_abstraction_score_laser:.1f}")
    print(f"Total relative D-removal score: {num_cd_bonds} (bonds) * {base_reactivity_cd:.1f} (reactivity) = {d_abstraction_score_laser:.1f}")

    # Calculate the probability of removing an H atom
    prob_h_removal_laser = h_abstraction_score_laser / total_reactivity_laser

    print(f"Result: The overall reaction is accelerated (total score from {total_reactivity_no_laser:.1f} to {total_reactivity_laser:.1f}).")
    print(f"The likelihood of H atom removal is now {prob_h_removal_laser:.1%}.\n")

    print("CONCLUSION: The simulation shows that exciting the C-H bond greatly accelerates the reaction and enhances the likelihood of H atom removal over D atoms.")

if __name__ == '__main__':
    simulate_reaction()