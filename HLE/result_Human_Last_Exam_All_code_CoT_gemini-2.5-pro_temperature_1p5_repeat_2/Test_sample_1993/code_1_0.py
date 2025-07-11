def diagnose_molecular_abnormality():
    """
    Analyzes patient data to determine the likely molecular abnormality.
    """

    # Patient's key clinical and genetic information
    age = 18
    height_percentile_range = "10th - 15th"
    karyotype = "Normal (46,XX)"
    symptoms = [
        "Amenorrhea and infertility",
        "Ovarian dysgenesis",
        "Short stature (persistently in the {} percentile)".format(height_percentile_range),
        "Fatigue and shortness of breath during physical activity",
        "Occasional elevated blood pressure"
    ]

    print("Step 1: Analyzing the patient's clinical presentation.")
    print("The patient is an {}-year-old female with the following key features:".format(age))
    for symptom in symptoms:
        print("- " + symptom)
    print("\nThis constellation of features (short stature, ovarian failure, and potential cardiac issues) is highly characteristic of Turner Syndrome.")

    print("\nStep 2: Considering the genetic study results.")
    print("The patient's karyotype was found to be: {}.".format(karyotype))
    print("This is a critical finding because it rules out classic Turner Syndrome, which is caused by a chromosomal abnormality (45,X0).")

    print("\nStep 3: Deducing the molecular abnormality.")
    print("Given the Turner-like phenotype but a normal chromosomal complement, the underlying cause is likely a single-gene mutation that phenocopies Turner Syndrome.")
    print("Noonan Syndrome is the most common genetic condition that presents with features overlapping Turner Syndrome, such as short stature and heart defects, and it can occur in individuals with a 46,XX karyotype.")
    print("The most frequent molecular cause of Noonan Syndrome (found in approximately 50% of cases) is a mutation in the PTPN11 gene, which is part of the RAS-MAPK signaling pathway critical for development.")

    print("\n--- Conclusion ---")
    print("The most likely molecular abnormality in this patient is a mutation in the PTPN11 gene.")


if __name__ == "__main__":
    diagnose_molecular_abnormality()
