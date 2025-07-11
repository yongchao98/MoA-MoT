import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak
    of orthorhombic NaMgH3.
    """
    # Step 1: Define crystal structure and lattice parameters
    print("This script calculates the Q-space position of the second major diffraction peak of NaMgH3.")
    print("\nStep 1: Crystal Information")
    print("NaMgH3 at room temperature has an orthorhombic crystal structure (space group Pnma).")
    a = 5.536  # in Angstroms
    b = 7.788  # in Angstroms
    c = 5.488  # in Angstroms
    print("Lattice parameters are:")
    print(f"  a = {a} Å")
    print(f"  b = {b} Å")
    print(f"  c = {c} Å")

    # Step 2: Identify the Miller indices (h,k,l) for the second major peak
    print("\nStep 2: Identify the Diffraction Peak")
    print("Based on diffraction data, the second major peak for this material corresponds to the (111) reflection.")
    h, k, l = 1, 1, 1
    print(f"The Miller indices are (h, k, l) = ({h}, {k}, {l}).")

    # Step 3: Calculate the d-spacing
    print("\nStep 3: Calculate the d-spacing")
    print("The formula for d-spacing in an orthorhombic system is: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
    
    # Calculate each term of the equation
    term_h = (h / a)**2
    term_k = (k / b)**2
    term_l = (l / c)**2
    
    # Calculate 1/d^2
    d_inv_sq = term_h + term_k + term_l
    print(f"  1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
    print(f"  1/d^2 = {term_h:.5f} + {term_k:.5f} + {term_l:.5f} = {d_inv_sq:.5f} Å⁻²")

    # Calculate d
    d = math.sqrt(1 / d_inv_sq)
    print(f"  d = 1 / sqrt({d_inv_sq:.5f}) = {d:.4f} Å")

    # Step 4: Calculate the Q-space position
    print("\nStep 4: Calculate the Q-space Position")
    print("The Q-space position is given by the formula: Q = 2 * pi / d")
    
    # Calculate Q
    q_value = 2 * math.pi / d
    print(f"  Q = (2 * {math.pi:.4f}) / {d:.4f}")
    
    print("\n--- FINAL RESULT ---")
    print(f"The second major diffraction peak ({h},{k},{l}) is located at Q = {q_value:.3f} Å⁻¹.")
    print("(Note: The measurement wavelength of 0.2952 Å is not needed for this calculation).")

    return q_value

# Run the calculation and store the final answer
final_q_value = calculate_q_space_position()
# The final answer is wrapped in <<<>>>
# print(f"\n<<<{final_q_value:.3f}>>>")
# The problem asks for the answer directly returned.

if __name__ == '__main__':
    # This block is for execution. The final response will have the value directly.
    pass

# Final answer to be returned by the AI model
# The code calculates it, let's put the value here.
# q_value = 1.803