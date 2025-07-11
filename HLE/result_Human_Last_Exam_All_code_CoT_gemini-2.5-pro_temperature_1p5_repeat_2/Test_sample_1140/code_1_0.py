import sys

def calculate_peptide_synthesis_yields():
    """
    Calculates and compares the theoretical yields for synthesizing a long peptide
    using Direct Solid-Phase Peptide Synthesis (SPPS) vs. Native Chemical Ligation (NCL).
    """

    # --- Parameters ---
    peptide_length = 100
    # Average coupling efficiency for each amino acid addition in SPPS.
    # 99.0% is a realistic efficiency for a well-optimized synthesis.
    spps_coupling_efficiency = 0.99
    
    # For NCL, we split the peptide into two fragments.
    fragment_1_length = 50
    fragment_2_length = 50 # Must sum to peptide_length
    
    # Efficiency of the ligation reaction itself. 70% is a reasonable estimate.
    ncl_ligation_efficiency = 0.70

    print("--- Comparing Peptide Synthesis Strategies ---")
    print(f"Target: A {peptide_length} amino acid peptide.")
    print(f"Assumption: A per-step coupling efficiency of {spps_coupling_efficiency:.1%} for SPPS.\n")

    # --- Calculation for Direct SPPS ---
    print("--- Method 1: Direct Solid-Phase Peptide Synthesis (SPPS) ---")
    direct_spps_steps = peptide_length - 1
    # The final yield is the efficiency multiplied by itself for each coupling step.
    direct_spps_yield = spps_coupling_efficiency ** direct_spps_steps
    
    print(f"This method requires {direct_spps_steps} sequential coupling steps.")
    print("The final theoretical yield is calculated as:")
    print(f"Yield = (Coupling Efficiency) ^ (Number of Steps)")
    print(f"      = {spps_coupling_efficiency} ^ {direct_spps_steps}")
    print(f"      = {direct_spps_yield:.3f}")
    print(f"\nResult: The theoretical yield for direct SPPS is {direct_spps_yield:.1%}.")
    print("The primary issue is the ~{:.1%} of impurities that are very difficult to remove.\n".format(100 - direct_spps_yield * 100))

    # --- Calculation for Native Chemical Ligation (NCL) ---
    print("--- Method 2: Native Chemical Ligation (NCL) ---")
    print(f"This method uses two shorter fragments ({fragment_1_length}aa and {fragment_2_length}aa) that are synthesized separately and then joined.")
    
    # Calculate yield for fragment 1
    frag1_spps_steps = fragment_1_length - 1
    frag1_spps_yield = spps_coupling_efficiency ** frag1_spps_steps
    print("\nFirst, calculate the SPPS yield for Fragment 1:")
    print(f"Yield = {spps_coupling_efficiency} ^ {frag1_spps_steps}")
    print(f"      = {frag1_spps_yield:.3f}  (or {frag1_spps_yield:.1%})")

    # The overall yield of the NCL process is the yield of the limiting purified fragment 
    # multiplied by the ligation efficiency. Assuming equal length fragments, we use one.
    overall_ncl_yield = frag1_spps_yield * ncl_ligation_efficiency
    
    print(f"\nNext, we ligate the purified fragments (ligation efficiency assumed to be {ncl_ligation_efficiency:.0%}).")
    print("The final theoretical yield is calculated as:")
    print(f"Overall Yield = (Fragment SPPS Yield) x (Ligation Efficiency)")
    print(f"              = {frag1_spps_yield:.3f} x {ncl_ligation_efficiency:.3f}")
    print(f"              = {overall_ncl_yield:.3f}")

    print(f"\nResult: The theoretical yield for the NCL strategy is {overall_ncl_yield:.1%}.")
    print("The key advantage is much higher purity due to easier purification of the shorter fragments.\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if overall_ncl_yield > direct_spps_yield:
        print(f"NCL is the superior method, offering a higher theoretical yield ({overall_ncl_yield:.1%})")
        print(f"compared to direct SPPS ({direct_spps_yield:.1%}) and resulting in a significantly purer product.")
    else:
        # This case is unlikely with these parameters but good to have
        print("Even with comparable theoretical yields, the vast improvement in product purity makes NCL the preferred method for long peptides.")

if __name__ == '__main__':
    calculate_peptide_synthesis_yields()