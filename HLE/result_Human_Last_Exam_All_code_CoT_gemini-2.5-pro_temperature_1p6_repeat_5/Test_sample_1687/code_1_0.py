def solve_medical_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis by evaluating
    key findings against the provided options.
    """

    print("Analyzing the clinical case based on key findings...\n")

    # --- Evidence Summary ---
    # 1. Procedure: "difficult" colonoscopy, "No polypectomy"
    # 2. Key Symptoms: "upper abdominal pain", "left-shoulder discomfort", "left upper quadrant" (LUQ) tenderness
    # 3. Signs of Hemorrhage: Hemoglobin drop from 11.7 to 6.5 g/dL, tachycardia, hypotension, conjunctival pallor
    # 4. Abdominal Exam: Worsening distension, guarding, peritoneal signs

    # --- Evaluating the Options ---

    # Option A: Colonic Perforation
    # Perforation can cause peritonitis and distension. However, the dominant
    # clinical picture is massive hemorrhage, not sepsis. The pain localization
    # to the LUQ and left shoulder is also less typical for a colonic perforation.
    perforation_score = 1  # For peritoneal signs
    print("Diagnosis A: Colonic perforation")
    print(" - Score: 1 (for peritoneal signs)")
    print(" - Rationale: While possible, the primary signs point to massive hemorrhage, not sepsis from perforation. Pain localization is atypical.\n")


    # Option B: Lower GI Bleeding
    # This involves bleeding from the colon. While the patient has massive bleeding,
    # the location of the pain (LUQ, left shoulder) is not characteristic.
    # A colonic bleed large enough to cause this shock would typically present with hematochezia (not mentioned).
    lower_gi_bleed_score = 1 - 1 # +1 for bleeding, -1 for atypical location
    print("Diagnosis B: Lower GI bleeding")
    print(" - Score: 1 (for massive hemorrhage) - 1 (for atypical pain location) = 0")
    print(" - Rationale: The pain in the LUQ and left shoulder strongly points away from a colonic source.\n")


    # Option C: Splenic Laceration
    # This is a known, though rare, complication of colonoscopy, especially a difficult one.
    # It perfectly explains the triad of:
    # 1. Left Upper Quadrant (LUQ) pain (location of the spleen).
    # 2. Left shoulder pain (Kehr's sign from diaphragmatic irritation by blood).
    # 3. Signs of severe hemorrhage (massive hemoglobin drop, shock).
    splenic_laceration_score = 1 + 1 + 1 # +1 for LUQ pain, +1 for left shoulder pain, +1 for severe hemorrhage
    print("Diagnosis C: Splenic laceration")
    print(" - Score: 1 (for LUQ pain) + 1 (for left shoulder pain) + 1 (for profound hemorrhage) = 3")
    print(" - Rationale: This diagnosis perfectly aligns with the location of pain, Kehr's sign, and the presentation of hemorrhagic shock after a difficult colonoscopy.\n")


    # Option D: Postpolypectomy Syndrome
    # This is directly ruled out by the case description.
    postpolypectomy_syndrome_score = 0 # Ruled out
    print("Diagnosis D: Postpolypectomy syndrome")
    print(" - Score: 0")
    print(" - Rationale: The text explicitly states 'No polypectomy was performed'.\n")

    # --- Final Conclusion ---
    print("="*40)
    print("Final Equation:")
    print(f"Score for (A) Colonic perforation = {perforation_score}")
    print(f"Score for (B) Lower GI bleeding = {lower_gi_bleed_score}")
    print(f"Score for (C) Splenic laceration = {splenic_laceration_score}")
    print(f"Score for (D) Postpolypectomy syndrome = {postpolypectomy_syndrome_score}")
    print("="*40)
    print("\nThe evidence most strongly supports Splenic laceration as the diagnosis.")

solve_medical_case()