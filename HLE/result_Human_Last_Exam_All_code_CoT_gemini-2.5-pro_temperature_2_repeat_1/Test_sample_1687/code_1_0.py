def diagnose_colonoscopy_complication():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    after a difficult colonoscopy.
    """
    # Key findings from the clinical case
    case_findings = {
        "procedure": "Difficult colonoscopy, no polypectomy",
        "pain_location": "Left upper quadrant (LUQ) and epigastrium",
        "referred_pain": "Left-shoulder discomfort (Kehr's sign)",
        "initial_hgb": 11.7, # g/dL
        "post_deterioration_hgb": 6.5, # g/dL
        "hemodynamics": "Developed tachycardia and hypotension (shock)",
        "exam_findings": "Worsening abdominal distension, guarding, peritoneal signs"
    }

    print("Analyzing the clinical case based on key findings...\n")

    # --- Reasoning Logic ---
    print("Evaluating Answer Choice A: Colonic Perforation")
    print(" - Possible, but the pain is specifically in the LUQ with referred left-shoulder pain. This specific pattern is less typical for a perforation.\n")

    print("Evaluating Answer Choice B: Lower GI Bleeding")
    print(" - This would explain the hemoglobin drop from {} g/dL to {} g/dL and shock.".format(case_findings['initial_hgb'], case_findings['post_deterioration_hgb']))
    print(" - However, the primary symptoms are severe LUQ pain and peritoneal signs, not overt bleeding from the rectum (hematochezia).\n")

    print("Evaluating Answer Choice D: Postpolypectomy Syndrome")
    print(" - This is definitively ruled out because the case explicitly states: '{}'.".format(case_findings['procedure']))
    print(" - Therefore, this diagnosis is impossible.\n")

    print("Evaluating Answer Choice C: Splenic Laceration")
    print(" - This diagnosis aligns perfectly with the classic clinical triad:")
    print("   1. Mechanism: A difficult colonoscopy, which can cause trauma to the nearby spleen.")
    print("   2. Location of Symptoms: Severe pain in the '{}', the anatomical location of the spleen.".format(case_findings['pain_location']))
    print("   3. Referred Pain: '{}', a classic sign (Kehr's sign) of diaphragmatic irritation from splenic hemorrhage.".format(case_findings['referred_pain']))
    print("   4. Hemorrhagic Shock: Explained by the acute drop in hemoglobin from {} to {} g/dL and subsequent tachycardia and hypotension.".format(case_findings['initial_hgb'], case_findings['post_deterioration_hgb']))
    print("   5. Abdominal Signs: '{}' are all consistent with a large volume of blood in the peritoneal cavity.".format(case_findings['exam_findings']))
    print("\nConclusion: Splenic Laceration provides the most comprehensive explanation for all the patient's signs and symptoms.")

    final_answer = "C"
    print("\nThe most likely diagnosis is Splenic Laceration.")
    print(f'<<<{final_answer}>>>')

diagnose_colonoscopy_complication()