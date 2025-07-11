def diagnose_patient():
    """
    Analyzes clinical findings to determine the most likely diagnosis.
    """
    # Step 1: Define the key clinical findings from the case.
    patient_findings = {
        "Age": "2-year-old",
        "Mass": "Pelvic mass (often abdominal/renal in origin)",
        "Vitals": "Hypertension (elevated blood pressure)",
        "Congenital Anomaly": "Aniridia (absence of the iris)",
        "Development": "Delayed speech",
        "Other": ["Conjunctival pallor (anemia)", "10th percentile weight/height (growth delay)"]
    }

    print("Analyzing patient presentation...")
    for key, value in patient_findings.items():
        print(f"- {key}: {value}")
    print("-" * 20)

    # Step 2: Identify the classic syndrome associated with these findings.
    print("Step 2: Identifying a potential syndrome.")
    print("The combination of Aniridia, a renal Mass, and Developmental Delay is highly characteristic of a specific syndrome.")
    print("This points towards WAGR Syndrome:")
    print("  W - Wilms Tumor (the most common type of Nephroblastoma)")
    print("  A - Aniridia")
    print("  G - Genitourinary anomalies")
    print("  R - Range of developmental delays")
    print("The patient exhibits at least three key features of WAGR syndrome.")
    print("-" * 20)

    # Step 3: Evaluate the differential diagnosis based on the syndrome.
    print("Step 3: Evaluating the answer choices.")
    print("A. Germ cell tumor: Not typically associated with aniridia or WAGR syndrome.")
    print("B. Astrocytoma: A brain tumor; does not present with a pelvic/abdominal mass.")
    print("C. Neuroblastoma: Can cause hypertension and an abdominal mass, but the link to aniridia is not a classic feature. It is a key differential but less likely given the aniridia.")
    print("D. Nephroblastoma (Wilms Tumor): This is a kidney tumor. It is the 'W' in WAGR syndrome. The patient's presentation of a mass, hypertension, and aniridia is a textbook case for Nephroblastoma associated with WAGR.")
    print("E. Ewing sarcoma: A bone/soft tissue tumor, not characteristically linked to this constellation of findings.")
    print("-" * 20)
    
    # Step 4: Conclude the most likely diagnosis.
    final_diagnosis = "D. Nephroblastoma"
    print("Conclusion: Based on the strong association with WAGR syndrome, the patient's pelvic mass is most likely a Nephroblastoma (Wilms tumor).")
    print(f"The correct answer is {final_diagnosis.split('.')[0]}")


if __name__ == "__main__":
    diagnose_patient()