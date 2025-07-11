def diagnose_sellar_mass(patient_age, patient_sex, tumor_location, cytologic_features):
    """
    This function simulates the diagnostic process for a sellar mass based on clinical and pathological findings.
    """
    
    print("Patient Clinical Information:")
    print(f"- Age: {patient_age} years")
    print(f"- Sex: {patient_sex}")
    print(f"- Tumor Location: {tumor_location}\n")
    
    print("Key Cytologic Findings from Crush Smear:")
    for feature in cytologic_features:
        print(f"- {feature}")
        
    # Diagnostic Logic
    diagnosis = "Undetermined"
    if (tumor_location == "sella turcica" and
        "Motononous cells" in cytologic_features and
        "'Salt-and-pepper' chromatin" in cytologic_features and
        "Eosinophilic granular cytoplasm" in cytologic_features):
        diagnosis = "Pituitary Adenoma"
        
    print("\n--- Diagnosis ---")
    print(f"The combination of a mass in the {tumor_location} in a {patient_age}-year-old patient")
    print("with cytologic features of uniform endocrine cells is characteristic of a Pituitary Adenoma.")
    
    return diagnosis

# Input data from the case
age = 28
sex = "Female"
location = "sella turcica"
findings = [
    "Cohesive sheets and single cells",
    "Motononous cells with round/oval nuclei",
    "'Salt-and-pepper' chromatin",
    "Eosinophilic granular cytoplasm",
    "Absence of significant pleomorphism or mitotic activity"
]

# Run the diagnostic simulation
final_diagnosis = diagnose_sellar_mass(age, sex, location, findings)

# The final answer is derived from this logical process.
# print(f"\nFinal Answer: {final_diagnosis}")