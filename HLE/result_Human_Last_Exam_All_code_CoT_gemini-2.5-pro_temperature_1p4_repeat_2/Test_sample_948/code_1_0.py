import numpy as np

def calculate_spin_correlation(state):
    """
    Calculates the nearest-neighbor spin correlation S_i^z * S_{i+1}^z.
    State representation: 1 for spin-up, -1 for spin-down, 0 for empty.
    S_z is proportional to these values, so we compute state[i] * state[i+1].
    """
    correlation = 0
    num_pairs = 0
    # Iterate through adjacent sites
    for i in range(len(state) - 1):
        # Only consider pairs of occupied sites
        if state[i] != 0 and state[i+1] != 0:
            correlation += (state[i]/2.0) * (state[i+1]/2.0)
            num_pairs += 1
    
    # Return average correlation per pair if pairs exist
    return correlation / num_pairs if num_pairs > 0 else 0

# --- State Examples ---

# A ferromagnetic-like state with spin segregation.
# This state is long-lived because opposite spins are far apart.
fm_state = [1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1]

# An anti-ferromagnetic-like state.
# This state is short-lived because opposite spins are adjacent,
# allowing for easy doublon formation and loss.
afm_state = [1, -1, 1, -1, 1, -1, 1, -1]


# Calculate correlations
fm_correlation = calculate_spin_correlation(fm_state)
afm_correlation = calculate_spin_correlation(afm_state)

print("Analysis of long-lived vs. short-lived states:")
print("-" * 50)
print(f"Long-lived (Ferromagnetic-like) State: {fm_state}")
print(f"Average nearest-neighbor spin correlation: {fm_correlation:.4f}")
print("This state has positive (ferromagnetic-like) spin correlations.\n")

print(f"Short-lived (Anti-ferromagnetic-like) State: {afm_state}")
print(f"Average nearest-neighbor spin correlation: {afm_correlation:.4f}")
print("This state has negative (anti-ferromagnetic-like) spin correlations.")
print("-" * 50)
print("\nConclusion:")
print("The loss mechanism eliminates states with anti-ferromagnetic character,")
print("leading to a system dominated by ferromagnetic-like correlations before it fully decays.")
print("\nTherefore, the properties of the system in the long-time limit are:")
print("1) Zero tunneling")
print("2) Zero particles")
print("3) Zero losses")
print("6) Ferromagnetic-like spin correlations")

p1, p2, p3, p6 = 1, 2, 3, 6
print("\nThe final equation from the property numbers is:")
print(f"{p1} + {p2} + {p3} + {p6} = {p1 + p2 + p3 + p6}")