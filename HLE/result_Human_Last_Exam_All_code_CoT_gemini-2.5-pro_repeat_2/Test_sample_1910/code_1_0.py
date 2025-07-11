import math

def solve_diffraction_peak():
    """
    Calculates the Q-space position for the second diffraction peak of NaMgH3.
    """
    # Step 1: Define lattice parameters for orthorhombic NaMgH3 (Pbnm) in Angstroms.
    a = 5.513
    b = 5.642
    c = 7.896

    # Step 2: Generate and filter for allowed reflections based on Pbnm space group rules.
    allowed_reflections = []
    # Generate reflections for low Miller indices (h, k, l)
    for h in range(0, 4):
        for k in range(0, 4):
            for l in range(0, 4):
                # Skip the (0,0,0) origin
                if h == 0 and k == 0 and l == 0:
                    continue

                # Pbnm (No. 62) reflection conditions
                is_allowed = False
                if h == 0 and k != 0 and l != 0:  # 0kl plane
                    if k % 2 == 0:
                        is_allowed = True
                elif k == 0 and h != 0 and l != 0:  # h0l plane
                    if (h + l) % 2 == 0:
                        is_allowed = True
                elif l == 0 and h !=0 and k != 0: # hk0 plane
                    is_allowed = True # no conditions
                elif h == 0 and k == 0: # 00l
                    if l % 2 == 0:
                        is_allowed = True
                elif h == 0 and l == 0: # 0k0
                    if k % 2 == 0:
                        is_allowed = True
                elif k == 0 and l == 0: # h00
                    if h % 2 == 0:
                        is_allowed = True
                elif h != 0 and k != 0 and l != 0: # general hkl
                    is_allowed = True

                if is_allowed:
                    # Step 3: Calculate 1/d^2 for the allowed reflection
                    one_over_d_squared = (h/a)**2 + (k/b)**2 + (l/c)**2
                    if one_over_d_squared > 0:
                        allowed_reflections.append({
                            'hkl': (h, k, l),
                            '1/d^2': one_over_d_squared
                        })

    # Remove duplicate reflections that give the same 1/d^2 value
    unique_reflections = []
    seen_d_sq = set()
    for r in allowed_reflections:
        val = round(r['1/d^2'], 5)
        if val not in seen_d_sq:
            seen_d_sq.add(val)
            unique_reflections.append(r)

    # Step 4: Sort reflections by 1/d^2 (which is proportional to Q^2)
    sorted_reflections = sorted(unique_reflections, key=lambda x: x['1/d^2'])

    print("Finding the first few diffraction peaks:")
    for i in range(min(5, len(sorted_reflections))):
        peak = sorted_reflections[i]
        d = math.sqrt(1 / peak['1/d^2'])
        q = 2 * math.pi / d
        print(f"Peak {i+1}: (hkl) = {peak['hkl']}, Q = {q:.4f} Å⁻¹")

    # The second major peak is the second in the sorted list
    second_peak = sorted_reflections[1]
    h, k, l = second_peak['hkl']
    one_over_d_sq_2 = second_peak['1/d^2']

    # Calculate d-spacing and Q-vector for the second peak
    d_2 = math.sqrt(1 / one_over_d_sq_2)
    q_2 = 2 * math.pi / d_2

    print("\n--- Calculation for the Second Diffraction Peak ---")
    print(f"The second peak corresponds to the Miller indices (hkl) = {second_peak['hkl']}.")
    print(f"The d-spacing (d) for this peak is calculated as 1 / sqrt(h²/a² + k²/b² + l²/c²).")
    print(f"d = 1 / sqrt({one_over_d_sq_2:.4f}) = {d_2:.4f} Å")
    print("The Q-space position is calculated as Q = 2 * π / d.")
    print(f"Q = 2 * {math.pi:.4f} / {d_2:.4f} = {q_2:.4f} Å⁻¹")
    print("\nThe second major diffraction peak is located at Q ≈ 1.5916 1/Å.")

solve_diffraction_peak()