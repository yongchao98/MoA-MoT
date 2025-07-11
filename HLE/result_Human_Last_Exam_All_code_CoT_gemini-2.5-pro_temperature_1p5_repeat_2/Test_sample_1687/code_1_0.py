def analyze_clinical_case():
    """
    Analyzes the provided clinical case to determine the most likely diagnosis.
    This function will print the step-by-step reasoning.
    """
    
    print("Step 1: Identify Key Findings from the Case")
    print("---------------------------------------------")
    print("Patient: 65-year-old female after a 'difficult' colonoscopy.")
    print("Procedure Note: 'No polypectomy was performed'.")
    print("Symptoms: Left upper quadrant (LUQ) and epigastric pain, plus left-shoulder discomfort (Kehr's sign).")
    print("Signs of Hemorrhagic Shock: Tachycardia (heart rate 105), hypotension (blood pressure 100/65), and dizziness.")
    print("Laboratory Evidence: Dramatic drop in hemoglobin from 11.7 g/dL to 6.5 g/dL, confirming massive blood loss.")
    print("Physical Exam: Worsening LUQ tenderness, abdominal distension, and peritoneal signs (guarding).\n")
    
    print("Step 2: Evaluate Each Diagnostic Option")
    print("------------------------------------------")
    print("A. Colonic Perforation: Possible, but less likely. Usually presents with signs of infection and is less likely to cause such rapid, massive hemorrhage. The specific pain pattern (LUQ + shoulder) is less characteristic.")
    print("B. Lower GI Bleeding: Unlikely. This diagnosis explains the blood loss but does not fit the pain pattern (LUQ + shoulder pain). It typically presents with rectal bleeding.")
    print("C. Splenic Laceration: Highly likely. This diagnosis perfectly explains all findings:")
    print("   - The spleen is in the LUQ, matching the pain location.")
    print("   - Bleeding from the spleen irritates the diaphragm, causing referred pain to the left shoulder (Kehr's sign).")
    print("   - The spleen is highly vascular, and an injury explains the massive hemorrhage (Hgb drop from 11.7 to 6.5) and shock.")
    print("   - It is a known, though rare, complication of difficult colonoscopies.")
    print("D. Postpolypectomy Syndrome: Ruled out. The case explicitly states 'No polypectomy was performed'.\n")

    print("Step 3: Conclusion")
    print("------------------")
    print("Based on the unique combination of LUQ pain, referred left-shoulder pain, and evidence of massive intra-abdominal hemorrhage after a difficult colonoscopy, splenic laceration is the most likely diagnosis.")
    print("\nFinal Answer Choice: C")

# Execute the analysis
analyze_clinical_case()