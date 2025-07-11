def solve_medical_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis.
    """

    # --- Step 1: Define Clinical Findings from the Case ---
    procedure_type = "Difficult screening colonoscopy"
    polypectomy_performed = False
    pain_location = "Left upper quadrant and epigastrium"
    referred_pain = "Left-shoulder discomfort"  # This is Kehr's sign.
    initial_hemoglobin = 11.7  # g/dL
    later_hemoglobin = 6.5    # g/dL
    post_transfusion_hemoglobin = 7.1 # g/dL
    hemodynamic_status = "Developed tachycardia and hypotension"
    physical_exam_findings = ["Worsening abdominal distension", "Guarding", "Peritoneal signs"]

    print("Analyzing the clinical case based on the following key findings:")
    print(f"- Procedure: {procedure_type}")
    print(f"- Polypectomy Performed: {polypectomy_performed}")
    print(f"- Pain Location: {pain_location}")
    print(f"- Referred Pain: {referred_pain}")
    print(f"- Hemodynamic Status: {hemodynamic_status}")
    print(f"- Hemoglobin Drop: From {initial_hemoglobin} g/dL to {later_hemoglobin} g/dL, confirming massive hemorrhage.\n")

    # --- Step 2: Evaluate Each Diagnosis ---

    print("Evaluating potential diagnoses:\n")

    # A. Colonic Perforation
    print("A. Colonic Perforation:")
    print("   - While possible after a difficult colonoscopy, this diagnosis is less likely.")
    print("   - The primary presentation is severe hemorrhage (hemoglobin drop from 11.7 to 6.5), not just peritonitis.")
    print("   - It does not typically explain the referred pain to the left shoulder (Kehr's sign).\n")

    # B. Lower GI Bleeding
    print("B. Lower GI Bleeding:")
    print("   - This is inconsistent with the findings.")
    print("   - Lower GI bleeding typically presents with blood per rectum (hematochezia), which was not reported.")
    print(f"   - The pain location ({pain_location}) is not typical for a colonic bleed.\n")

    # C. Splenic Laceration
    print("C. Splenic Laceration:")
    print("   - This diagnosis is highly consistent with all major findings.")
    print(f"   - Mechanism: A known, though rare, complication of a '{procedure_type}'.")
    print(f"   - Symptoms: Perfectly explains the '{pain_location}' and '{referred_pain}'.")
    print(f"   - Signs of Hemorrhage: Explains the '{hemodynamic_status}' and the dramatic hemoglobin drop from {initial_hemoglobin} to {later_hemoglobin}.")
    print(f"   - Abdominal Signs: Consistent with the '{physical_exam_findings}' due to blood accumulation (hemoperitoneum).\n")
    
    # D. Postpolypectomy Syndrome
    print("D. Postpolypectomy Syndrome:")
    print("   - This diagnosis can be ruled out completely.")
    print(f"   - The case explicitly states that a polypectomy was performed: {polypectomy_performed}.\n")

    # --- Step 3: Conclude ---
    print("Conclusion: The combination of a difficult colonoscopy, LUQ pain, Kehr's sign (left shoulder pain), and profound hemorrhagic shock makes Splenic Laceration the most likely diagnosis.")

solve_medical_case()
print("<<<C>>>")