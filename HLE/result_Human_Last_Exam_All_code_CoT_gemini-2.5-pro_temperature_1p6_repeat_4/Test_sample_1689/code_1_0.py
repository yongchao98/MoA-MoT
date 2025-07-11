import textwrap

def solve_clinical_case():
    """
    Analyzes the clinical case and determines the best next step in diagnosis.
    """
    
    # Patient presentation summary
    patient_age = 43
    symptom = "Itchy rash in both axillae"
    key_history = "Started wearing new workout clothes 3 weeks prior"
    key_finding = "Rash on the periphery of the axillary folds, sparing the vaults"
    provisional_diagnosis = "Allergic contact dermatitis from textiles (clothing)"

    # The question to be answered
    question = "Which of the following is the best next step in diagnosis?"
    
    # Answer choices
    choices = {
        'A': 'Skin biopsy',
        'B': 'KOH preparation',
        'C': 'Topical steroid',
        'D': 'Patch test',
        'E': 'None of the above'
    }

    # Logical analysis
    analysis = f"""
1. The patient's history and physical exam provide strong clues. The {patient_age}-year-old man developed a rash after starting a new workout regimen involving new clothes.

2. The location of the rash is critical. The rash is on the posterior axillary folds, where friction from clothing is highest, and it spares the axillary vault. This pattern is characteristic of textile (clothing) contact dermatitis, not deodorant dermatitis (which would primarily affect the vault).

3. The question asks for the next best DIAGNOSTIC step, not treatment.
   - Choice A (Skin biopsy) is invasive and not the first-line test for suspected contact dermatitis.
   - Choice B (KOH preparation) is used to diagnose fungal infections, which are less likely given the specific history and rash distribution.
   - Choice C (Topical steroid) is a TREATMENT, not a diagnostic test.
   - Choice D (Patch test) is the gold-standard diagnostic tool to identify the specific allergen causing allergic contact dermatitis. It is the most appropriate step to confirm the suspected diagnosis of textile dermatitis.

4. Therefore, the patch test is the best next step to confirm the diagnosis.
"""

    correct_answer_key = 'D'
    correct_answer_text = choices[correct_answer_key]

    print("### Clinical Case Analysis ###")
    print(textwrap.fill(f"The case describes a {patient_age}-year-old man with an itchy rash. The key information pointing to the diagnosis is the recent history of new workout clothes and the specific rash distribution, which spares the axillary vault.", 80))
    print("\n### Evaluating the Options ###")
    print(textwrap.fill(analysis, 80))
    print("-" * 30)
    print(f"Conclusion: The best next step in diagnosis is D: {correct_answer_text}.")
    print("-" * 30)

solve_clinical_case()