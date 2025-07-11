import pandas as pd

def analyze_lc_data():
    """
    Analyzes the provided liquid crystal data to answer the two-part question.
    """
    # Step 1: Analyze the effect on relaxation dynamics by comparing data from the plots.
    # Data extracted for ring 1 (blue diamonds) from both plots.
    # Temperatures are in Kelvin (K), and relaxation times <τ> are in nanoseconds (ns).
    n1_ring1_data = {
        325: 200,
        340: 50,
        350: 25,
        375: 8,
        400: 3
    }
    m1_ring1_data = {
        325: 150,
        340: 30,
        350: 15,
        375: 9,
        400: 3
    }

    print("--- Part 1 Analysis: Effect on Relaxation Dynamics ---")
    print("Comparing relaxation time <τ> for ring 1 in N1 (nonmethylated) vs. M1 (methylated):")
    
    faster_count = 0
    slower_count = 0
    
    for temp, tau_n1 in n1_ring1_data.items():
        tau_m1 = m1_ring1_data.get(temp)
        if tau_m1 is not None:
            if tau_m1 < tau_n1:
                comparison = "faster"
                faster_count += 1
            else:
                comparison = "slower or same"
                slower_count += 1
            print(f"At {temp} K: N1 τ = {tau_n1} ns, M1 τ = {tau_m1} ns. Methylated ring dynamics are {comparison}.")

    print("\nConclusion for Part 1:")
    if faster_count > slower_count:
        print("The data shows that the relaxation time <τ> for the methylated ring is generally lower.")
        print("This means the addition of a methyl group leads to a DECREASED relaxation time (faster dynamics).")
    else:
        print("The data shows that the relaxation time <τ> for the methylated ring is generally higher.")
        print("This means the addition of a methyl group leads to an INCREASED relaxation time (slower dynamics).")

    # Step 2: Analyze the effect on the nematic-isotropic transition temperature (T_NI).
    print("\n--- Part 2 Analysis: Effect on Nematic-Isotropic Transition Temperature (T_NI) ---")
    print("Physical Chemistry Principle: Adding a bulky lateral group (like a methyl group) to a liquid crystal molecule disrupts the ordered packing.")
    print("This disruption weakens the intermolecular forces that stabilize the nematic phase.")
    print("\nConclusion for Part 2:")
    print("A less stable nematic phase requires less thermal energy to become a disordered isotropic liquid.")
    print("Therefore, the addition of a methyl group is expected to DECREASE the nematic-isotropic transition temperature.")

    # Step 3: Combine conclusions to select the final answer.
    print("\n--- Final Answer Selection ---")
    print("Summary of conclusions:")
    print("1. Relaxation time: DECREASED.")
    print("2. Nematic-isotropic transition temperature: DECREASED.")
    print("\nMatching this with the answer choices leads to option E.")
    
    answer_e = "E. 1. Because of steric bulk, the addition of a methyl group to the structure decreases the correlation time of the methylated ring  relative to the nonmethylated mesogen, leading to a decreased relaxation time. 2. The addition of a methyl group disrupts crystallinity, leading to a lower nematic-isotropic transition temperature."
    print(f"\nFull text of the correct answer:\n{answer_e}")

# Execute the analysis
analyze_lc_data()

print("\n<<<E>>>")