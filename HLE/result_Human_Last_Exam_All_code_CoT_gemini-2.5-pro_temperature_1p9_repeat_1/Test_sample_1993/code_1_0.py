def diagnose_patient_condition():
    """
    Analyzes a patient's clinical profile to determine the likely molecular abnormality.
    This script simulates the diagnostic reasoning process for the given medical case.
    """
    # Patient's key clinical and genetic findings
    patient_findings = {
        "Karyotype": "46,XX (Normal Female)",
        "Primary Complaint": "Amenorrhea and Infertility",
        "Ultrasound Finding": "Ovarian Dysgenesis",
        "Other Symptoms": ["Short Stature", "Exercise Intolerance", "Intermittent Hypertension"]
    }

    print("--- Diagnostic Reasoning Process ---")

    # Step 1: Initial Assessment
    print("\nStep 1: Evaluating the core clinical picture.")
    print(f"The patient presents with features resembling Turner Syndrome, such as {patient_findings['Other Symptoms'][0]} and {patient_findings['Ultrasound Finding']}.")
    print("However, the crucial piece of information is the karyotype.")

    # Step 2: Interpreting the Karyotype
    print("\nStep 2: Analyzing the genetic report.")
    print(f"The patient's karyotype is '{patient_findings['Karyotype']}'.")
    if patient_findings["Karyotype"] == "46,XX (Normal Female)":
        print("This rules out classic Turner Syndrome, which is caused by a 45,X karyotype.")
        print("The diagnosis is therefore a condition that mimics Turner Syndrome but occurs in individuals with a normal chromosomal complement.")
    else:
        print("The provided karyotype does not match the case description.")
        return

    # Step 3: Identifying the specific condition
    print("\nStep 3: Naming the condition.")
    print(f"The combination of a normal female karyotype and '{patient_findings['Ultrasound Finding']}' points to '46,XX Gonadal Dysgenesis'.")
    print("The next step is to find the underlying molecular cause.")

    # Step 4: Pinpointing the Molecular Abnormality
    print("\nStep 4: Determining the likely molecular abnormality.")
    print("We are looking for a gene mutation that can cause this condition.")
    print("A key candidate is the BMP15 (Bone Morphogenetic Protein 15) gene because:")
    print("  1. It is located on the X-chromosome, which can explain the overlap with some Turner Syndrome features.")
    print("  2. It plays a critical role in ovarian follicle development.")
    print("  3. Mutations in BMP15 are a known cause of ovarian dysgenesis, primary ovarian insufficiency, and infertility in 46,XX females.")

    # Final Conclusion
    final_diagnosis = "A mutation in the BMP15 gene"
    print("\n--- Final Conclusion ---")
    print(f"The patient's presentation of amenorrhea, infertility, and ovarian dysgenesis, despite a normal 46,XX karyotype, makes the most likely molecular abnormality:")
    print(f">>> {final_diagnosis} <<<")

# Run the diagnostic script
diagnose_patient_condition()