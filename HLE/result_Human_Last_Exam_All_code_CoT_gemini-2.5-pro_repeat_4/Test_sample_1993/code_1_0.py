def analyze_patient_case():
    """
    Analyzes the clinical vignette of a patient with Turner-like features
    but a normal karyotype, and identifies the likely molecular abnormality.
    """

    # Patient data from the prompt
    patient_age = 18
    height_percentile_lower = 10
    height_percentile_upper = 15
    chromosome_complement = "46,XX"

    print("--- Patient Case Analysis ---")
    print(f"Patient Profile: An {patient_age}-year-old female.")
    print(f"Key Findings:")
    print(f"  - Stature: Persistently in the {height_percentile_lower}th - {height_percentile_upper}th percentile.")
    print(f"  - Symptoms: Amenorrhea, infertility, fatigue, and elevated blood pressure.")
    print(f"  - Physical Exam: Ovarian dysgenesis observed on ultrasound.")
    print(f"  - Genetics: Reported normal chromosomal complement of {chromosome_complement}.")
    print("\n--- Diagnostic Reasoning ---")
    print("1. The combination of short stature, ovarian dysgenesis, and potential cardiovascular issues is classic for Turner Syndrome.")
    print(f"2. However, the patient's reported normal '{chromosome_complement}' karyotype rules out classic Turner Syndrome (which is 45,X).")
    print("3. We must consider a molecular cause that explains all symptoms while being consistent with a 46-chromosome count.")
    print("4. The features of Turner Syndrome are primarily caused by having only one functional copy of critical genes on the short arm of the X chromosome (Xp).")
    print("5. A molecular abnormality that would cause this, yet might be missed on a standard karyotype, is a small deletion of genetic material.")
    print("\n--- Conclusion ---")
    print("The most likely molecular abnormality is a deletion on the short arm of one of the X chromosomes [del(Xp)].")
    print("This single abnormality explains:")
    print("  - Short Stature: Due to the loss of one copy of the SHOX gene.")
    print("  - Ovarian Dysgenesis: Due to the loss of other critical genes for ovarian function on Xp.")
    print("  - Turner-like Phenotype: The overall presentation is accounted for by this loss of genetic material.")
    print("  - 'Normal' Karyotype: A small deletion can be missed by standard chromosome analysis, leading to a report of 46 chromosomes.")

if __name__ == "__main__":
    analyze_patient_case()
