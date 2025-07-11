def diagnose_sellar_mass(clinical_data, histo_findings):
    """
    Analyzes clinical and histological data to provide a diagnosis.

    Args:
        clinical_data (dict): A dictionary with patient's clinical information.
        histo_findings (dict): A dictionary with key histological features.

    Returns:
        str: The final diagnosis.
    """
    # Key features for Mucormycosis diagnosis
    mucor_features = {
        "hyphae_width": "broad",
        "septation": "pauciseptate",
        "branching": "right-angle"
    }

    # Compare observed findings with classic features of Mucormycosis
    if (histo_findings.get("organism_type") == "fungal hyphae" and
        histo_findings.get("hyphae_width") == mucor_features["hyphae_width"] and
        histo_findings.get("septation") == mucor_features["septation"] and
        histo_findings.get("branching") == mucor_features["branching"]):
        
        diagnosis = "Mucormycosis"
    else:
        diagnosis = "Further investigation needed"
        
    # Print the analysis
    print("Clinical and Pathological Analysis:")
    print("-" * 35)
    print(f"Patient: {clinical_data['age']} year old {clinical_data['gender']}")
    print(f"Location: Mass in the {clinical_data['location']}")
    print("\nHistological Findings:")
    print(f"- Organism Type: {histo_findings['organism_type']}")
    print(f"- Hyphae Width: {histo_findings['hyphae_width']}")
    print(f"- Septation: {histo_findings['septation']}")
    print(f"- Branching Pattern: {histo_findings['branching']}")
    print("-" * 35)
    print(f"\nConclusion: The findings are characteristic of {diagnosis}.")
    
    return diagnosis

# --- Input Data ---
# Clinical information provided
patient_info = {
    "age": 28,
    "gender": "female",
    "location": "sella turcica"
}

# Observations from the provided histopathology image
histology_observations = {
    "organism_type": "fungal hyphae",
    "hyphae_width": "broad",
    "septation": "pauciseptate",
    "branching": "right-angle"
}

# Run the diagnostic function
final_diagnosis = diagnose_sellar_mass(patient_info, histology_observations)