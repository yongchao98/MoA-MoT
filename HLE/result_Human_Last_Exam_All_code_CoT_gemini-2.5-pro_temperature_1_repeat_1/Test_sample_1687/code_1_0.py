def diagnose_patient():
    """
    Analyzes the clinical findings to determine the most likely diagnosis.
    """
    print("Analyzing the clinical case based on the provided findings...")
    print("-" * 50)

    # Key findings from the case study
    procedure = "Difficult colonoscopy, no polypectomy"
    pain_location = "Left upper quadrant (LUQ) and left-shoulder discomfort"
    initial_hemoglobin = 11.7
    lowest_hemoglobin = 6.5
    hemodynamic_status = "Developed tachycardia and hypotension (shock)"

    print(f"Patient History: {procedure}")
    print(f"Key Symptoms: {pain_location}")
    print(f"Evidence of Bleeding: Hemoglobin dropped from {initial_hemoglobin} g/dL to {lowest_hemoglobin} g/dL.")
    print(f"Hemodynamic Status: {hemodynamic_status}")

    print("\nEvaluating the potential diagnoses:")
    print("A. Colonic perforation: Possible, but the specific combination of LUQ pain and referred left-shoulder pain is less typical. Hemorrhage is usually not this severe unless a major vessel is involved.")
    print("B. Lower GI bleeding: This typically presents with blood in the stool and is less likely to cause severe LUQ and shoulder pain.")
    print("D. Postpolypectomy syndrome: This is ruled out because no polypectomy was performed.")
    
    print("\nAnalyzing the most likely option:")
    print("C. Splenic laceration:")
    print("  - It is a known complication of difficult colonoscopies due to traction on the splenocolic ligament.")
    print("  - It directly explains the LUQ pain (location of the spleen).")
    print("  - Bleeding into the abdomen irritates the diaphragm and phrenic nerve, causing referred pain to the left shoulder (Kehr's sign). This is a classic sign.")
    print(f"  - The spleen is a highly vascular organ, and an injury explains the rapid and severe blood loss, as shown by the hemoglobin drop from {initial_hemoglobin} to {lowest_hemoglobin}, leading to shock.")

    print("\n" + "-" * 50)
    print("Conclusion: The combination of a difficult colonoscopy, LUQ pain, left-shoulder pain (Kehr's sign), and profound hemorrhagic shock makes Splenic Laceration the most compelling diagnosis.")

diagnose_patient()