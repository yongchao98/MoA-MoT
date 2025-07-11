def solve_medical_case():
    """
    This function analyzes the clinical case to determine the most likely diagnosis.
    """
    # Key clinical findings from the text
    procedure = "Difficult screening colonoscopy"
    polypectomy_performed = False
    pain_location_1 = "Upper abdominal pain and left-shoulder discomfort"
    pain_location_2 = "Tenderness to palpation involving her left upper quadrant and epigastrium"
    initial_hemoglobin = 11.7
    post_deterioration_hemoglobin = 6.5
    hemodynamic_status = "Developed tachycardia and hypotension (hemorrhagic shock)"
    peritoneal_signs = "Developed a significant amount of guarding in these areas... demonstrated peritoneal signs"

    print("Analyzing the clinical findings to determine the most likely diagnosis:")
    print("-" * 60)
    print(f"1. Procedure Context: The patient had a '{procedure}'. This can lead to complications like perforation or injury to adjacent organs.")
    print(f"2. Pain Localization: The pain is consistently localized to the 'left upper quadrant' and radiates to the 'left shoulder'. This specific shoulder pain is known as Kehr's sign and strongly suggests irritation of the diaphragm, often by blood.")
    print(f"3. Evidence of Bleeding: The hemoglobin dropped sharply from {initial_hemoglobin} g/dL to {post_deterioration_hemoglobin} g/dL, and the patient became hemodynamically unstable. This confirms significant internal hemorrhage.")
    print("-" * 60)
    print("Evaluating the answer choices:")
    print("A. Colonic perforation: Unlikely. While it causes abdominal pain, it does not typically cause referred left shoulder pain. The pain localization points away from this diagnosis.")
    print("B. Lower GI bleeding: Unlikely. This would explain blood loss but not the upper abdominal and left shoulder pain.")
    print("D. Postpolypectomy syndrome: Impossible. The case explicitly states that no polypectomy was performed.")
    print("C. Splenic laceration: Most Likely. The spleen is in the left upper quadrant. Injury is a known, though rare, complication of difficult colonoscopies. This diagnosis perfectly explains all key symptoms:")
    print("   - Pain in the left upper quadrant (location of the spleen).")
    print("   - Referred left shoulder pain (Kehr's sign from sub-diaphragmatic blood).")
    print("   - Massive internal hemorrhage and shock from a ruptured vascular organ.")
    print("-" * 60)
    print("Conclusion: The combination of a difficult colonoscopy, LUQ pain, Kehr's sign, and hemorrhagic shock makes splenic laceration the most likely diagnosis.")

solve_medical_case()
print("<<<C>>>")