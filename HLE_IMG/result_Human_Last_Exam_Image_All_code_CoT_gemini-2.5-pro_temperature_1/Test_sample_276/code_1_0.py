def solve_medical_case():
    """
    This function analyzes the provided medical case and determines the most likely diagnosis.
    """
    # Patient Data
    age = 67
    gender = "female"
    symptoms = ["low-grade fever", "weight loss", "fatigue", "diarrhea", "acute abdominal pain"]
    duration_months = 1
    history = ["gallstones", "uveitis", "arthritis", "animal shelter volunteer"]
    exam_findings = ["right hemiabdomen tenderness", "guarding"]
    lab_wbc = 13000
    lab_fobt = "positive"
    ct_findings = "marked ileocecal wall thickening"

    # Analysis
    # The key features are subacute onset (1 month) of constitutional and GI symptoms,
    # localization to the right lower quadrant (ileocecal area), history of uveitis/arthritis,
    # and a specific exposure risk (animal shelter).
    # The CT confirms ileocolitis.

    # Evaluating top choices:
    # A. Crohn's Disease: Fits symptoms and location, but age of onset is less typical.
    # C. Ileocecal Tuberculosis: A classic mimic, but no specific TB risk factors are mentioned.
    # B. Yersinia Colitis: Fits symptoms, location, and CT findings. Crucially, it directly
    #    links to the specific risk factor of animal exposure. Yersinia can cause a
    #    subacute presentation and is known to trigger reactive arthritis and uveitis.

    # Conclusion: The most likely diagnosis is the one that ties all the information together,
    # especially the unique risk factor.
    most_likely_diagnosis_code = 'B'
    most_likely_diagnosis_name = "Yersinia Colitis"

    print(f"Patient Presentation Summary:")
    print(f"- Age: {age}")
    print(f"- Symptoms: {', '.join(symptoms)}")
    print(f"- Duration: {duration_months} month")
    print(f"- Key History: {', '.join(history)}")
    print(f"- Exam: {', '.join(exam_findings)}")
    print(f"- WBC Count: {lab_wbc}")
    print(f"- CT Findings: {ct_findings}")
    print("\nAnalysis:")
    print("The constellation of symptoms, subacute course, ileocecal involvement on CT,")
    print("and history of inflammatory conditions like uveitis and arthritis points towards")
    print("an inflammatory or infectious process in the right lower quadrant.")
    print("The specific history of volunteering at an animal shelter is a strong clue for a")
    print("zoonotic infection. Yersinia enterocolitica is a bacterium transmitted from animals")
    print("that classically causes ileocolitis, mimicking Crohn's disease, and can be associated")
    print("with reactive arthritis and uveitis.")
    print("\nConclusion:")
    print(f"The most likely diagnosis is ({most_likely_diagnosis_code}) {most_likely_diagnosis_name}.")

solve_medical_case()