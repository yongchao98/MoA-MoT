import math

def calculate_q_space():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # Step 1: Define crystal structure parameters and Miller indices
    # Lattice parameters for orthorhombic NaMgH3 at room temperature (in Angstrom)
    # Source: Schouwink et al., J. Phys. Chem. C 2012, 116, 10829–10836
    a = 5.4957
    b = 7.8447
    c = 5.4800

    # Miller indices (h, k, l) for the second major diffraction peak, which is the (121) reflection.
    h, k, l = 1, 2, 1

    print("This script calculates the Q-space position of the second major diffraction peak of NaMgH3.")
    print(f"The crystal structure is orthorhombic with lattice parameters: a={a} Å, b={b} Å, c={c} Å.")
    print(f"The second major peak corresponds to the ({h},{k},{l}) reflection.\n")

    # Step 2: Calculate the interplanar spacing 'd'
    # The formula for an orthorhombic system is: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
    d_inv_sq = (h / a)**2 + (k / b)**2 + (l / c)**2
    d = math.sqrt(1 / d_inv_sq)

    # Step 3: Calculate the Q-space value
    # The formula is Q = 2 * pi / d
    pi = math.pi
    Q = (2 * pi) / d

    # Step 4: Print the detailed calculation and the final result
    print("--- Calculation Steps ---")
    print("1. Calculate the interplanar spacing (d):")
    print(f"   1/d² = ({h}/{a})² + ({k}/{b})² + ({l}/{c})² = {d_inv_sq:.5f} Å⁻²")
    print(f"   d = 1 / sqrt({d_inv_sq:.5f}) = {d:.4f} Å")
    
    print("\n2. Calculate the Q-space position (Q):")
    print("   Q = 2 * π / d")
    # Outputting each number in the final equation as requested
    print(f"   Q = 2 * {pi:.5f} / {d:.4f} = {Q:.3f} Å⁻¹")
    
    print("\n--- Final Answer ---")
    print(f"The second major diffraction peak of NaMgH3 is located at Q = {Q:.3f} 1/Angstrom.")

if __name__ == '__main__':
    calculate_q_space()