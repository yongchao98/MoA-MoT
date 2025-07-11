def find_diagnosis():
    """
    This function analyzes the patient's symptoms to determine the most likely diagnosis
    by identifying a classic syndrome.
    """
    # Patient's key clinical features
    features = {
        "Age": "2-year-old",
        "Mass": "Pelvic mass",
        "Eye Finding": "Aniridia",
        "Development": "Delayed speech",
        "Vitals": "Elevated blood pressure (hypertension)"
    }

    print("Analyzing the combination of patient's symptoms:")
    for key, value in features.items():
        print(f"- {key}: {value}")

    print("\n--- Diagnostic Reasoning ---")
    print("The presence of Aniridia alongside a pelvic mass and developmental delay in a young child is a classic presentation.")
    print("This constellation of findings strongly suggests WAGR Syndrome, which is an acronym for:")
    print("W - Wilms Tumor (also known as Nephroblastoma)")
    print("A - Aniridia")
    print("G - Genitourinary anomalies")
    print("R - Retardation (mental/developmental)")

    print("\nThe patient's pelvic mass is consistent with a Wilms Tumor (Nephroblastoma).")
    print("The patient's hypertension is also a common feature of Wilms Tumor, often caused by the tumor secreting renin.")

    print("\n--- Conclusion ---")
    print("Comparing this clinical picture to the answer choices, Nephroblastoma is the most fitting diagnosis.")
    print("The correct choice is D.")


# Run the diagnostic analysis
find_diagnosis()