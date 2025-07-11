import numpy as np

def calculate_frictional_force(normal_load, alignment_factor):
    """
    Calculates the frictional force in a simplified model of atomic friction,
    demonstrating the principles of superlubricity.

    Friction is dependent on the normal load and the alignment of atomic lattices.
    
    Args:
        normal_load (float): The force pressing the surfaces together (in nanoNewtons).
        alignment_factor (float): A value from 0.0 (incommensurate, perfect superlubricity)
                                  to 1.0 (commensurate, high friction).
    """
    # Base friction coefficient for an ideal superlubric state (close to zero)
    mu_base = 0.001
    
    # A factor representing how much friction increases due to lattice alignment (coherence)
    mu_coherence_effect = 0.3
    
    # Effective friction coefficient increases with alignment
    effective_mu = mu_base + (alignment_factor * mu_coherence_effect)
    
    # Frictional force is the effective coefficient times the normal load
    frictional_force = effective_mu * normal_load
    
    print(f"Calculating friction with Normal Load = {normal_load} nN and Alignment Factor = {alignment_factor}:")
    print(f"Effective Coefficient (μ) = Base μ [{mu_base}] + (Alignment Factor [{alignment_factor}] * Coherence Effect [{mu_coherence_effect}])")
    print(f"Effective Coefficient (μ) = {mu_base} + {alignment_factor * mu_coherence_effect:.3f} = {effective_mu:.3f}")
    print(f"Frictional Force = Effective μ [{effective_mu:.3f}] * Normal Load [{normal_load}]")
    print(f"Resulting Frictional Force = {frictional_force:.3f} nN\n")

# --- Demonstration ---

# Case 1: Low alignment (superlubric state)
print("--- Case 1: Superlubric State ---")
load_1 = 10.0  # nN
alignment_1 = 0.05 # Almost incommensurate
calculate_frictional_force(load_1, alignment_1)

# Case 2: High alignment (high friction state)
print("--- Case 2: High Friction State ---")
load_2 = 10.0 # Same load
alignment_2 = 0.95 # Almost commensurate
calculate_frictional_force(load_2, alignment_2)

# Case 3: Increased normal load in the high friction state
print("--- Case 3: High Friction with Higher Load ---")
load_3 = 50.0 # Increased load
alignment_3 = 0.95 # Same high alignment
calculate_frictional_force(load_3, alignment_3)
