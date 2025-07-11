def solve_case():
    """
    This function analyzes the clinical case and determines the most likely diagnosis.
    """
    # Patient Profile
    age = 67
    sex = "female"
    
    # Key Symptoms (present over 1 month)
    symptoms = ["low-grade fever", "weight loss", "fatigue", "diarrhea"]
    acute_symptom = "acute abdominal pain"
    
    # Past Medical History
    history = ["gallstones", "uveitis", "arthritis"]
    
    # Exam Findings
    exam = ["right hemiabdomen tenderness", "guarding"]
    
    # Lab Results
    wbc_count = 13000  # cells/mcL
    fecal_occult_blood = "positive"
    
    # Imaging Findings (CT Scan)
    ct_findings = ["marked wall thickening of the ileocecal region"]
    
    # Analysis
    # The constellation of findings is highly suggestive of a specific diagnosis.
    # 1. Location: Right lower quadrant pain + ileocecal thickening on CT points to disease of the terminal ileum and/or cecum.
    # 2. Chronicity: 1-month history of symptoms suggests a chronic inflammatory process, not an acute infection or event.
    # 3. Extraintestinal Manifestations: The history of uveitis and arthritis are classic extraintestinal manifestations of Inflammatory Bowel Disease (IBD).
    # 4. Synthesis: Crohn's Disease is the form of IBD that classically affects the ileocecal region ("ileocolitis") and is strongly associated with both uveitis and arthritis. The patient's age is consistent with a late-onset presentation.
    
    diagnosis = "A. Crohn's Disease"
    
    print("Patient Profile:")
    print(f" - Age: {age}")
    print(f" - Sex: {sex}")
    
    print("\nKey Clinical Features:")
    print(f" - Chronic Symptoms: {', '.join(symptoms)}")
    print(f" - Past Medical History: {', '.join(history)}")
    print(f" - Exam Findings: {', '.join(exam)}")
    print(f" - Lab Findings: WBC {wbc_count}, FOBT {fecal_occult_blood}")
    print(f" - CT Findings: {ct_findings[0]}")
    
    print("\nConclusion:")
    print("The combination of chronic ileocecal inflammation with classic extraintestinal manifestations (uveitis, arthritis) makes Crohn's Disease the most likely diagnosis.")
    
    # Final Answer format is not code-related, but for clarity:
    # print("\nFinal Answer Choice: A")

solve_case()