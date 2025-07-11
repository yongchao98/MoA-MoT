def analyze_phage_experiments():
    """
    Analyzes the results of two experiments to determine the roles of
    the bacterial RP system and the phage's XY operon.
    """
    # Data from Experiment 1: Plaque-Forming Units (CFU) per microliter
    # In bacteria without the RP defense system
    cfu_wt_no_rp = 100000
    cfu_deltaXY_no_rp = 100000

    # In bacteria with the RP defense system
    cfu_wt_with_rp = 80000
    cfu_deltaXY_with_rp = 40000

    print("--- Analysis of Experiment 1: Phage Virulence ---")

    # Part 1: Quantify the resistance provided by the RP system.
    # Resistance is the reduction in phage efficiency (CFU).
    # Calculate for the wild-type phage (PhageDE3-wt)
    resistance_vs_wt = ((cfu_wt_no_rp - cfu_wt_with_rp) / cfu_wt_no_rp) * 100
    print(f"Resistance effect of RP system against PhageDE3-wt:")
    print(f"({cfu_wt_no_rp} - {cfu_wt_with_rp}) / {cfu_wt_no_rp} * 100 = {resistance_vs_wt:.0f}% reduction in CFU.")
    
    # Calculate for the mutant phage (PhageDE3-deltaXY)
    resistance_vs_deltaXY = ((cfu_deltaXY_no_rp - cfu_deltaXY_with_rp) / cfu_deltaXY_no_rp) * 100
    print(f"\nResistance effect of RP system against PhageDE3-deltaXY:")
    print(f"({cfu_deltaXY_no_rp} - {cfu_deltaXY_with_rp}) / {cfu_deltaXY_no_rp} * 100 = {resistance_vs_deltaXY:.0f}% reduction in CFU.")
    print("\nConclusion from Exp 1: The RP system increases bacterial resistance against the phage.")

    print("\n--- Analysis of Experiment 2: Mass Spectrometry ---")
    print("The 500 Da molecule is detected ONLY when:")
    print("1. Phage has the XY operon (PhageDE3-wt)")
    print("2. Bacterium has the RP system")
    print("\nConclusion from Exp 2: The product of the XY enzymes is only synthesized in the presence of the RP system.")

    print("\n--- Overall Conclusion ---")
    print("The results from both experiments show that the RP system provides resistance, and the phage's counter-defense (operon XY) is specifically activated by the presence of this RP system to produce its active molecule. This matches statement H.")

analyze_phage_experiments()
<<<H>>>