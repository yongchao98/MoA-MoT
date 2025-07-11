def analyze_clinical_case():
    """
    Analyzes the patient's clinical data to determine the likely molecular abnormality.
    """
    # Patient Clinical Data
    patient_age = 18
    height_percentile_lower = 10
    height_percentile_upper = 15
    karyotype = "46,XX (Normal)"
    
    symptoms = {
        "Reproductive": "Amenorrhea, infertility, ovarian dysgenesis",
        "Growth": f"Short stature (persistently in {height_percentile_lower}th-{height_percentile_upper}th percentile)",
        "Cardiopulmonary": "Fatigue and shortness of breath on exertion, elevated blood pressure"
    }
    
    # Print the analysis logic
    print("Step 1: Evaluating the Phenotype")
    print(f"The patient, age {patient_age}, presents with a classic Turner Syndrome phenotype:")
    for system, desc in symptoms.items():
        print(f"- {system}: {desc}")
    
    print("\nStep 2: Evaluating the Genotype")
    print(f"A standard karyotype analysis returned a normal result: {karyotype}.")
    print("This rules out classic monosomy X (45,X) as the cause.")

    print("\nStep 3: Synthesizing the Findings")
    print("The presence of a Turner Syndrome phenotype with a normal karyotype points to a submicroscopic genetic defect.")
    print("The genes responsible for stature and ovarian function in this syndrome are located on the short arm of the X chromosome (Xp).")
    
    print("\nConclusion:")
    print("A deletion in the Xp region would lead to haploinsufficiency of these critical genes, causing the observed symptoms without being detected by standard karyotyping.")
    
    final_diagnosis = "Deletion in the short arm of the X chromosome (Xp)"
    print(f"\nLikely Molecular Abnormality: {final_diagnosis}")

# Run the analysis
analyze_clinical_case()