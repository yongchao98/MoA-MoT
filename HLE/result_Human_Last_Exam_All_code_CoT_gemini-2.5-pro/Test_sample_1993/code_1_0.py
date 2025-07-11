def diagnose_patient():
    """
    Analyzes a clinical vignette to determine the likely molecular abnormality.
    """
    # Patient data from the clinical vignette
    patient_age = 18
    height_percentile_lower = 10
    height_percentile_upper = 15

    # Step 1: Summarize the clinical presentation
    print("--- Patient Presentation Summary ---")
    print(f"Age: {patient_age} years")
    print(f"Key Symptoms: Amenorrhea, infertility, fatigue, and shortness of breath.")
    print(f"Key Findings: Short stature ({height_percentile_lower}th-{height_percentile_upper}th percentile), ovarian dysgenesis, occasional high blood pressure, and a normal 46,XX karyotype.")
    print("-" * 34 + "\n")

    # Step 2: Print the diagnostic reasoning
    print("--- Diagnostic Analysis ---")
    print("1. The combination of short stature, ovarian dysgenesis, and cardiovascular signs strongly suggests Turner Syndrome.")
    print("2. However, the patient's normal 46,XX karyotype rules out classic Turner Syndrome (45,X).")
    print("3. We must therefore consider a 'phenocopy' - a condition that mimics Turner Syndrome but has a different underlying cause.")
    print("4. Noonan Syndrome is the most common condition with Turner-like features (short stature, heart defects) and a normal karyotype.")
    print("5. Although ovarian dysgenesis is an atypical feature of Noonan Syndrome, it is reported in some cases. It provides the most likely single genetic cause for the patient's entire set of symptoms.")
    print("-" * 27 + "\n")

    # Step 3: State the conclusion
    print("--- Conclusion ---")
    print("The likely diagnosis is Noonan Syndrome.")
    print("The most common molecular abnormality causing Noonan Syndrome is a mutation in the PTPN11 gene.")
    print("-" * 18)

# Execute the diagnostic function
diagnose_patient()