def provide_diagnosis():
    """
    This function provides the diagnosis based on the clinical and pathological findings.
    """
    patient_age = 28
    location = "sella turcica"
    
    # Key Cytological Findings:
    # 1. Monotonous population of cells
    # 2. Round to oval nuclei
    # 3. "Salt-and-pepper" stippled chromatin
    # 4. Syncytial arrangement (indistinct cell borders)
    # 5. Eosinophilic cytoplasm
    
    diagnosis = "Pituitary Adenoma"
    
    print(f"Patient Age: {patient_age}")
    print(f"Location of Mass: {location}")
    print("\nBased on the classic neuroendocrine cytological features (monotonous cells with 'salt-and-pepper' chromatin) in the context of a sellar mass, the diagnosis is:")
    print(f"\nDiagnosis: {diagnosis}")

provide_diagnosis()