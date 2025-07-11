def analyze_experiment_limitations():
    """
    This function explains the fundamental limitation of bulk calorimetric experiments
    for studying nucleic acid thermodynamics.
    """
    print("Analyzing the options for the limitation of bulk calorimetric experiments:")
    print("="*70)

    # Explanation for choices A and B
    print("Option A (ΔCp = 0) and B (T-independent parameters):")
    print("These are simplifying assumptions of the most common Nearest-Neighbor model, not inherent limitations of the experimental technique itself.")
    print("The experimental data from a bulk measurement could be fit to a more complex model where these parameters are not constant.")
    print("\n")

    # Explanation for choice C
    print("Option C (Impossibility to capture heterogeneity):")
    print("This is the key limitation of any 'bulk' experiment. A bulk measurement averages the signal from millions or billions of molecules.")
    print("In reality, the population is heterogeneous: not all molecules melt at the exact same time or follow the exact same pathway.")
    print("A bulk experiment cannot distinguish these individual molecular behaviors. It measures a population average, masking the underlying diversity and any rare intermediate states.")
    print("This is a fundamental consequence of ensemble averaging.")
    print("\n")
    
    # Explanation for choices D and E
    print("Option D (Temperature oscillations) and E (Temperature cannot be controlled):")
    print("These are factually incorrect. Modern calorimetry provides precise and stable temperature control and ramping, which is essential for the experiment to work.")
    print("="*70)

    print("Conclusion: The most fundamental limitation is that bulk experiments average out molecular differences, thus cannot capture heterogeneity.")

analyze_experiment_limitations()