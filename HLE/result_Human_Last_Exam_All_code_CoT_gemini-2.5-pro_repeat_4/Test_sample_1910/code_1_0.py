import math

def solve_diffraction_peak():
    """
    Calculates the Q-space position of the second major diffraction peak for NaMgH3.
    """
    # Step 1: Define crystal structure properties
    # Orthorhombic Pnma structure for NaMgH3 at room temperature
    a = 5.517  # Angstrom
    b = 7.778  # Angstrom
    c = 5.421  # Angstrom

    print("Analyzing the diffraction pattern for NaMgH3 (orthorhombic, Pnma).")
    print(f"Lattice parameters: a={a} Å, b={b} Å, c={c} Å")

    def check_reflection_allowed(h, k, l):
        """
        Checks if a reflection (h,k,l) is allowed for the Pnma space group.
        Conditions for Pnma (No. 62):
        - 0kl: k+l = 2n
        - h0l: h = 2n
        - hk0: h = 2n
        - h00: h = 2n
        - 0k0: k = 2n
        - 00l: l = 2n
        """
        if h == 0 and k == 0 and l == 0:
            return False
        
        # Check special conditions based on which indices are zero
        if h == 0 and k == 0:  # 00l
            return l % 2 == 0
        if h == 0 and l == 0:  # 0k0
            return k % 2 == 0
        if k == 0 and l == 0:  # h00
            return h % 2 == 0
        if h == 0:  # 0kl
            return (k + l) % 2 == 0
        if k == 0:  # h0l
            return h % 2 == 0
        if l == 0:  # hk0
            return h % 2 == 0
        
        # If none of the above, it's a general (h,k,l) reflection with no conditions
        return True

    # Step 2: Find all allowed peaks and their Q values
    allowed_peaks = []
    # Iterate through a reasonable range of Miller indices
    for h in range(0, 4):
        for k in range(0, 4):
            for l in range(0, 4):
                if check_reflection_allowed(h, k, l):
                    # Calculate Q for the allowed peak
                    d_sq_inv = (h / a)**2 + (k / b)**2 + (l / c)**2
                    if d_sq_inv > 0:
                        d_hkl = 1 / math.sqrt(d_sq_inv)
                        q_hkl = 2 * math.pi / d_hkl
                        allowed_peaks.append({'hkl': (h, k, l), 'Q': q_hkl})

    # Step 3: Sort the peaks by Q value
    sorted_peaks = sorted(allowed_peaks, key=lambda x: x['Q'])

    # Step 4: Identify the second peak and show the calculation
    if len(sorted_peaks) < 2:
        print("Could not find at least two diffraction peaks in the calculated range.")
        return

    second_peak = sorted_peaks[1]
    h, k, l = second_peak['hkl']
    
    print(f"\nThe first major peak is {sorted_peaks[0]['hkl']} at Q={sorted_peaks[0]['Q']:.4f} Å⁻¹.")
    print(f"The second major peak is {second_peak['hkl']} at Q={second_peak['Q']:.4f} Å⁻¹.")
    
    print("\n--- Detailed Calculation for the Second Major Peak ---")
    
    # Equation Part 1: 1/d^2
    print("\n1. Calculate the inverse square of the d-spacing (1/d²):")
    print(f"   Formula: 1/d² = (h/a)² + (k/b)² + (l/c)²")
    print(f"   Values:  1/d² = ({h}/{a})² + ({k}/{b})² + ({l}/{c})²")
    d_sq_inv_h = (h / a)**2
    d_sq_inv_k = (k / b)**2
    d_sq_inv_l = (l / c)**2
    d_sq_inv_total = d_sq_inv_h + d_sq_inv_k + d_sq_inv_l
    print(f"   Result:  1/d² = {d_sq_inv_h:.5f} + {d_sq_inv_k:.5f} + {d_sq_inv_l:.5f} = {d_sq_inv_total:.5f} Å⁻²")

    # Equation Part 2: d
    d_hkl = 1 / math.sqrt(d_sq_inv_total)
    print(f"\n2. Calculate the d-spacing (d):")
    print(f"   Formula: d = 1 / sqrt(1/d²)")
    print(f"   Values:  d = 1 / sqrt({d_sq_inv_total:.5f})")
    print(f"   Result:  d = {d_hkl:.4f} Å")
    
    # Equation Part 3: Q
    pi_val = math.pi
    q_final = 2 * pi_val / d_hkl
    print(f"\n3. Calculate the Q-space position (Q):")
    print(f"   Formula: Q = 2 * π / d")
    print(f"   Values:  Q = 2 * {pi_val:.5f} / {d_hkl:.4f}")
    print(f"   Result:  Q = {q_final:.4f} Å⁻¹")

    print(f"\nConclusion: The second major diffraction peak is located in Q-space at {q_final:.4f} 1/Å.")
    
    # Return the final answer in the specified format
    return q_final

if __name__ == '__main__':
    q_value = solve_diffraction_peak()
    # The final answer is directly returned after the main output.
    # print(f'<<<{q_value:.4f}>>>') # This line is for the final answer format
    
solve_diffraction_peak()