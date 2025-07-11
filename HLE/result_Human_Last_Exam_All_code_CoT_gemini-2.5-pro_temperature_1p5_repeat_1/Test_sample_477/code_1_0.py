def model_ctg_instability():
    """
    This function models the impact of LIG1 on CTG repeat instability.
    It calculates a hypothetical "instability score" based on whether LIG1 is active.
    """

    # --- Model Parameters ---
    # Baseline instability from repeat sequence alone
    baseline_instability = 5.0
    
    # A multiplier representing the expansion-prone repair pathway that requires LIG1.
    # A value > 1 means LIG1's pathway actively expands the repeats.
    lig1_dependent_expansion_factor = 3.0
    
    # A fixed, low level of instability from an alternative pathway when LIG1 is absent.
    alternative_pathway_instability = 1.0

    # --- Case 1: Normal Cell with active LIG1 ---
    print("--- Simulating Instability with Active LIG1 ---")
    
    # In this case, the main expansion-prone pathway is active.
    final_instability_with_lig1 = baseline_instability * lig1_dependent_expansion_factor
    
    print("Equation with LIG1:")
    print(f"Final Instability = Baseline Instability ({baseline_instability}) * LIG1-Dependent Expansion Factor ({lig1_dependent_expansion_factor})")
    print(f"Result: {baseline_instability} * {lig1_dependent_expansion_factor} = {final_instability_with_lig1}")
    print("\n")
    
    # --- Case 2: Cell with LIG1 Knockout ---
    print("--- Simulating Instability with LIG1 Knockout ---")
    
    # With LIG1 knocked out, the expansion-prone pathway is disabled.
    # Instability is now governed by the much less severe alternative pathway.
    final_instability_without_lig1 = alternative_pathway_instability

    print("Equation without LIG1:")
    print(f"Final Instability = Alternative Pathway Instability ({alternative_pathway_instability})")
    print(f"Result: {final_instability_without_lig1}")
    print("\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if final_instability_with_lig1 > final_instability_without_lig1:
        print(f"The instability score is reduced from {final_instability_with_lig1} to {final_instability_without_lig1} upon LIG1 knockout.")
        print("This demonstrates that knocking out LIG1 leads to reduced somatic instability.")
    else:
        print("The model does not show a reduction in instability.")

# Execute the model
model_ctg_instability()