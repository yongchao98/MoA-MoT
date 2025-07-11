def solve_oled_radical_problem():
    """
    This function explains the main disadvantage of air-stable organic radicals in OLEDs
    and determines the correct answer from the given choices.
    """
    explanation_title = "Analysis of Radicals in OLEDs"
    print(explanation_title)
    print("=" * len(explanation_title))

    # Step 1: Explain the advantage (context)
    print("\n[Step 1] The Advantage of Organic Radicals:")
    print("In OLEDs, electrically generated excitons are 25% singlets (light-emitting) and 75% triplets (dark).")
    print("Organic radicals have a doublet ground state, which allows them to harvest energy from both singlets and triplets.")
    print("This bypasses the 25% limit, theoretically allowing for 100% internal quantum efficiency (IQE).\n")

    # Step 2: Explain the primary disadvantage
    print("[Step 2] The Core Disadvantage: Quenching")
    print("The unpaired electron in a radical, while beneficial for harvesting excitons, also makes it an effective quencher.")
    print("An exciton (the excited state that should produce light) can interact with a ground-state radical molecule nearby.")
    print("This interaction provides a fast, non-radiative pathway for the exciton to lose its energy without producing light.")
    print("This process is called 'doublet-exciton quenching' and is a major efficiency loss mechanism.\n")

    # Step 3: Evaluate the choices based on the analysis
    print("[Step 3] Evaluating the Answer Choices:")
    print("A: Incorrect. The question specifies 'air-stable' radicals.")
    print("B: Incorrect. Wide emission spectrum is a color purity issue, not the fundamental efficiency problem.")
    print("C: Partially correct but is a consequence, not the root cause. Low luminance results from low efficiency.")
    print("D: Correct. This accurately identifies the root cause. Quenching of excitons by radicals directly lowers the External Quantum Efficiency (EQE), which is the primary performance metric.")
    print("E: Incorrect. Delocalization is a stabilizing strategy and not the direct cause of low EQE.\n")

    # Step 4: Final Conclusion
    print("[Conclusion]")
    print("The primary disadvantage is that the radical centers themselves act as quenchers for excitons, which severely lowers the device's External Quantum Efficiency (EQE).")

solve_oled_radical_problem()
<<<D>>>