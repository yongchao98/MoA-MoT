import math

def solve():
    """
    Calculates the Q-space position of the second major diffraction peak
    for NaMgH3 perovskite.
    """
    # Lattice parameters for orthorhombic NaMgH3 (Pnma) in Angstroms
    a = 5.534
    b = 7.784
    c = 5.488

    # The problem mentions a wavelength of 0.2952 Å.
    # Note: The Q-space position (Q = 2*pi/d) is independent of the wavelength.
    # Wavelength is only needed to calculate the scattering angle 2-theta.

    def is_allowed_pnma(h, k, l):
        """
        Checks if a reflection (h,k,l) is allowed for space group Pnma (#62).
        Systematic absences are:
        - 0kl: k+l must be even
        - h0l: h must be even
        """
        # Exclude the (0,0,0) point, which is not a reflection
        if h == 0 and k == 0 and l == 0:
            return False
        
        # Condition for 0kl reflections
        if h == 0 and (k + l) % 2 != 0:
            return False
            
        # Condition for h0l reflections
        if k == 0 and h % 2 != 0:
            return False
            
        return True

    # Dictionary to store unique Q-values and their corresponding Miller indices
    q_to_hkl = {}
    
    # Generate reflections for a range of Miller indices
    max_index = 4
    for h in range(max_index + 1):
        for k in range(max_index + 1):
            for l in range(max_index + 1):
                # We only need to consider one octant due to symmetry
                if is_allowed_pnma(h, k, l):
                    # Calculate 1/d^2
                    one_over_d_sq = (h / a)**2 + (k / b)**2 + (l / c)**2
                    
                    if one_over_d_sq > 1e-9: # Avoid division by zero or log(0)
                        q_val = 2 * math.pi * math.sqrt(one_over_d_sq)
                        
                        # Store the reflection, avoiding near-duplicates by rounding Q
                        q_rounded = round(q_val, 5)
                        if q_rounded not in q_to_hkl:
                            q_to_hkl[q_rounded] = (h, k, l)
    
    # Sort the unique Q-values in ascending order
    unique_q_values = sorted(q_to_hkl.keys())
    
    # --- Identify the second peak ---
    # The first peak may be a cluster of very close reflections.
    # We define a peak as a set of reflections whose Q-values are very close.
    # We find the first Q-value that is distinctly separate from the initial cluster.
    q_peak1 = unique_q_values[0]
    
    second_peak_q = 0
    second_peak_hkl = None
    
    # A tolerance of 1% of the first peak's Q value to distinguish peaks
    tolerance = 0.01 * q_peak1 
    
    for q in unique_q_values:
        if (q - q_peak1) > tolerance:
            second_peak_q = q
            second_peak_hkl = q_to_hkl[q]
            break

    # --- Print the detailed calculation for the second peak ---
    h, k, l = second_peak_hkl
    one_over_d_sq_val = (h/a)**2 + (k/b)**2 + (l/c)**2

    print(f"The second major diffraction peak corresponds to the Miller indices (h, k, l) = {second_peak_hkl}.")
    print("The Q-space position is calculated using the formula for an orthorhombic crystal:")
    print("\nQ = 2 * pi * sqrt((h/a)^2 + (k/b)^2 + (l/c)^2)")
    print(f"Q = 2 * {math.pi:.5f} * sqrt(({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2)")
    print(f"Q = 2 * {math.pi:.5f} * sqrt({(h/a)**2:.5f} + {(k/b)**2:.5f} + {(l/c)**2:.5f})")
    print(f"Q = 2 * {math.pi:.5f} * sqrt({one_over_d_sq_val:.5f})")
    print(f"Q = {second_peak_q:.4f} 1/Å")


solve()
<<<1.6151>>>