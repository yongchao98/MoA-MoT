import math

def calculate_diffraction_peaks():
    """
    Calculates and identifies the second major diffraction peak for NaMgH3.
    """
    # Step 1: Define crystal parameters for orthorhombic NaMgH3 (Pnma)
    a = 5.642  # Angstrom
    b = 7.915  # Angstrom
    c = 5.568  # Angstrom

    print("Crystal structure: Orthorhombic NaMgH3 (Pnma)")
    print(f"Lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å\n")

    # Step 2: Function to check Pnma reflection conditions
    def is_allowed_pnma(h, k, l):
        # (0,0,0) is not a reflection
        if h == 0 and k == 0 and l == 0:
            return False
        # For space group Pnma (No. 62):
        # 0kl reflection: k+l must be even
        if h == 0 and (k + l) % 2 != 0:
            return False
        # hk0 reflection: h must be even
        if l == 0 and h % 2 != 0:
            return False
        # h0l reflection: h+l must be even
        if k == 0 and (h + l) % 2 != 0:
            return False
        return True

    # Step 3: Generate and calculate Q-space positions for allowed reflections
    allowed_peaks = []
    # Search a reasonable range for Miller indices to find the first few peaks
    for h in range(5):
        for k in range(5):
            for l in range(5):
                if is_allowed_pnma(h, k, l):
                    # Formula for d-spacing in an orthorhombic system
                    d_inv_sq = (h / a)**2 + (k / b)**2 + (l / c)**2
                    if d_inv_sq > 0:
                        d = 1 / math.sqrt(d_inv_sq)
                        # Formula for Q-space position
                        Q = 2 * math.pi / d
                        allowed_peaks.append({'Q': Q, 'hkl': (h, k, l), 'd': d})

    # Sort peaks by increasing Q value
    allowed_peaks.sort(key=lambda p: p['Q'])

    # Step 4: Group nearby peaks and identify the second major peak group
    if not allowed_peaks:
        print("No peaks found.")
        return

    grouped_peaks = []
    current_group = [allowed_peaks[0]]
    for i in range(1, len(allowed_peaks)):
        # Group peaks with a Q difference less than 0.1 Å⁻¹
        if allowed_peaks[i]['Q'] - current_group[-1]['Q'] < 0.1:
            current_group.append(allowed_peaks[i])
        else:
            grouped_peaks.append(current_group)
            current_group = [allowed_peaks[i]]
    grouped_peaks.append(current_group)

    # The first few groups usually correspond to the "major" peaks
    first_peak_group = grouped_peaks[0]
    second_peak_group = grouped_peaks[1]

    # --- Outputting Results ---
    print("The first major peak group consists of the following reflections:")
    for peak in first_peak_group:
        print(f"  hkl: {peak['hkl']}, d-spacing: {peak['d']:.4f} Å, Q: {peak['Q']:.4f} 1/Å")

    print("\nThe second major peak group consists of the following reflections:")
    for peak in second_peak_group:
        print(f"  hkl: {peak['hkl']}, d-spacing: {peak['d']:.4f} Å, Q: {peak['Q']:.4f} 1/Å")
    
    # Calculate the average Q value for the second peak group
    q_sum_second_peak = sum(p['Q'] for p in second_peak_group)
    q_avg_second_peak = q_sum_second_peak / len(second_peak_group)
    
    # Show the calculation for a representative peak from the second group
    rep_peak = second_peak_group[0]
    h, k, l = rep_peak['hkl']
    d_rep = rep_peak['d']
    Q_rep = rep_peak['Q']

    print("\n--- Detailed Calculation for a Representative Peak in the Second Group ---")
    print(f"Let's use the first reflection in this group, hkl = {h,k,l}:")
    print("The d-spacing is calculated as:")
    print("d_hkl = 1 / sqrt((h/a)² + (k/b)² + (l/c)²)")
    print(f"d_{h}{k}{l} = 1 / sqrt(({h}/{a})² + ({k}/{b})² + ({l}/{c})²)")
    print(f"d_{h}{k}{l} = {d_rep:.4f} Å")
    print("\nThe Q-space position is calculated as:")
    print("Q = 2 * π / d")
    print(f"Q = 2 * {math.pi:.4f} / {d_rep:.4f}")
    print(f"Q = {Q_rep:.4f} 1/Å")
    
    print(f"\nThe average Q-space position for the second major diffraction peak is {q_avg_second_peak:.3f} 1/Å.")
    
    # Final answer in the required format
    print(f"\n<<<{q_avg_second_peak:.3f}>>>")

# Run the calculation
calculate_diffraction_peaks()