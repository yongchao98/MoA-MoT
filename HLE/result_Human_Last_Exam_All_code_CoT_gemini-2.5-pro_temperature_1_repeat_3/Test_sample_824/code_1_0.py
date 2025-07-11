def find_lab_parameter():
    """
    Analyzes the clinical scenario to identify the most relevant lab parameter.
    """
    patient_history = {
        "Long-term Symptoms": "Facial rash, joint pain, fever, hematuria (suggesting SLE)",
        "New Symptoms": "Dry eyes, dry mouth (suggesting secondary Sjögren's)",
        "Clinical Event": "Rebound inflammatory flare after stopping corticosteroids",
        "Outcome": "Rapid progression to end-stage kidney disease"
    }

    # The rapid renal decline in an SLE patient points to severe lupus nephritis.
    # Lupus nephritis is an immune-complex-mediated disease.
    # Immune complexes deposit in the kidneys and activate the complement system.
    # This activation consumes complement proteins.
    cause = "Active lupus nephritis driven by immune complex deposition."
    
    # Key complement proteins consumed in this process are C3 and C4.
    key_protein_1 = "C3"
    key_protein_2 = "C4"

    print("Based on the clinical scenario, the patient likely has Systemic Lupus Erythematosus (SLE).")
    print("The rapid deterioration of kidney function is characteristic of a severe flare of lupus nephritis.")
    print("\nThis condition is caused by immune complexes depositing in the kidneys, which triggers a strong inflammatory response.")
    print("This response actively consumes complement proteins from the blood.")
    print("\nTherefore, the lab parameter that best indicates the cause of the rapid renal decline would be low serum complement levels.")
    print(f"Specifically, measuring the levels of {key_protein_1} and {key_protein_2} would show a significant decrease, confirming the diagnosis.")

find_lab_parameter()