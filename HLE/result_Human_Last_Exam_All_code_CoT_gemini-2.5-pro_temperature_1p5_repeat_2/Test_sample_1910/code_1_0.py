import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define known parameters
    # Lattice parameter 'a' for cubic NaMgH3 in Angstroms (Å)
    a = 3.90
    # Miller indices (h, k, l) for the second major diffraction peak in a cubic perovskite
    h, k, l = 2, 0, 0

    print(f"The calculation is for NaMgH3, a cubic perovskite with lattice parameter a = {a} Å.")
    print(f"The second major diffraction peak corresponds to the (h, k, l) = ({h}, {k}, {l}) reflection.")
    print("-" * 30)

    # Step 2: Calculate the sum of the squares of the Miller indices
    hkl_sq_sum = h**2 + k**2 + l**2
    print(f"The sum of squares of Miller indices is: h² + k² + l² = {hkl_sq_sum}")

    # Step 3: Calculate the interplanar spacing (d-spacing)
    # Formula for a cubic system: d = a / sqrt(h² + k² + l²)
    d_hkl = a / math.sqrt(hkl_sq_sum)
    print(f"The interplanar spacing d_hkl is calculated as: d = {a} / sqrt({hkl_sq_sum}) = {d_hkl:.4f} Å")

    # Step 4: Calculate the Q-space position
    # Formula: Q = 2 * pi / d
    Q_hkl = (2 * math.pi) / d_hkl
    
    print("\nThe Q-space position is calculated as:")
    # Using 'eval' is generally unsafe, but here we construct the string ourselves for display purposes.
    # We are showing the numbers used in the final equation.
    final_equation = f"Q = (2 * {math.pi:.4f}) / {d_hkl:.4f} = {Q_hkl:.4f} Å⁻¹"
    print(final_equation)

if __name__ == "__main__":
    calculate_q_space_position()
