import numpy as np

def check_quantum_evolution(diagram_label, sz, s_plus_mag, s):
    """
    Checks if the given values from a quantum evolution plot are physically valid for a single qubit.

    Args:
        diagram_label (str): The label of the diagram (e.g., 'A').
        sz (float): The expectation value <σz>.
        s_plus_mag (float): The magnitude of the expectation value |<σ+>|.
        s (float): The von Neumann entropy S.
    """
    print(f"--- Checking Diagram {diagram_label} ---")
    valid = True

    # 1. Check <σz> bounds
    if not (-1 <= sz <= 1):
        print(f"Violation: <σz> = {sz} is outside the allowed range [-1, 1].")
        valid = False

    # 2. Check Entropy bounds
    max_entropy = np.log(2)
    if s < 0:
        print(f"Violation: Entropy S = {s} cannot be negative.")
        valid = False
    if s > max_entropy:
        print(f"Violation: Entropy S = {s:.3f} is greater than the maximum for a qubit, log(2) ≈ {max_entropy:.3f}.")
        valid = False

    # 3. Check Bloch sphere constraint
    bloch_constraint_value = 4 * (s_plus_mag**2) + sz**2
    if bloch_constraint_value > 1:
        # Re-print the equation with numbers
        print(f"Violation: Bloch sphere constraint failed.")
        print(f"4 * |<σ+>|^2 + <σz>^2 = 4 * {s_plus_mag}^2 + {sz}^2 = {bloch_constraint_value:.3f}, which is > 1.")
        valid = False
    
    if valid:
        print("This point appears to be physically valid.")
        # Print the check for the valid case as well
        print(f"Bloch sphere constraint check: 4 * |<σ+>|^2 + <σz>^2 = 4 * {s_plus_mag}^2 + {sz}^2 = {bloch_constraint_value:.3f}, which is <= 1.")

    print("-" * (20 + len(diagram_label)))


# Picking representative points from each plot to test for violations
# Diagram A: at t ≈ 1.5
check_quantum_evolution('A', sz=-0.2, s_plus_mag=0.9, s=0.3)

# Diagram B: at t ≈ 1
check_quantum_evolution('B', sz=0.65, s_plus_mag=0.65, s=0.1)

# Diagram C: at t ≈ 3
check_quantum_evolution('C', sz=1.7, s_plus_mag=0.6, s=-1.2)

# Diagram D: at t ≈ 6
check_quantum_evolution('D', sz=0.4, s_plus_mag=0.4, s=0.8)

# Diagram E: at t ≈ 1.5
check_quantum_evolution('E', sz=0.7, s_plus_mag=0.7, s=0.15)

# Diagram F: Checking two points to be thorough
# Point 1 (t ≈ 1): <σz> is max
check_quantum_evolution('F', sz=0.7, s_plus_mag=0.2, s=0.1)
# Point 2 (t ≈ 2.5): |<σ+>| is max
check_quantum_evolution('F', sz=0.5, s_plus_mag=0.3, s=0.2)
