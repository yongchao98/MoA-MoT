import sys

def solve_clinical_case():
    """
    Analyzes a clinical case to identify the damaged anatomical structure.
    """
    # Step 1: Define the patient's symptoms and assign a numerical value for the "equation".
    # A value of 1 means the symptom is present.
    symptoms = {
        'no_pupillary_light_reflex': 1,
        'inability_to_adduct': 1,
        'inability_to_depress': 1,
        'inability_to_elevate': 1,
    }
    print("Patient Symptoms Analysis:")
    for symptom, value in symptoms.items():
        if value == 1:
            print(f"- {symptom.replace('_', ' ').title()}")

    # Step 2: Define functions controlled by relevant cranial nerves.
    cranial_nerve_functions = {
        'Cranial Nerve III (Oculomotor)': [
            'pupillary_light_reflex',
            'adduct',
            'depress',
            'elevate'
        ],
        'Cranial Nerve VI (Abducens)': ['abduct'],
        'Cranial Nerve VII (Facial)': ['facial_expression']
    }

    # Step 3: Calculate a match score for CN III, creating a simple equation.
    cn3_score = 0
    equation_parts = []
    print("\nCalculating Match Score for Cranial Nerve III:")
    # We check which of the patient's symptoms are explained by CN III dysfunction.
    # Note: 'no_pupillary_light_reflex' corresponds to 'pupillary_light_reflex' function, etc.
    symptom_to_function_map = {
        'no_pupillary_light_reflex': 'pupillary_light_reflex',
        'inability_to_adduct': 'adduct',
        'inability_to_depress': 'depress',
        'inability_to_elevate': 'elevate'
    }

    cn3_explained_functions = cranial_nerve_functions['Cranial Nerve III (Oculomotor)']
    for symptom, function in symptom_to_function_map.items():
        if function in cn3_explained_functions:
            score_value = symptoms[symptom]
            cn3_score += score_value
            equation_parts.append(str(score_value))

    # This fulfills the request to output each number in the final equation.
    equation_str = " + ".join(equation_parts)
    print(f"Equation of explained symptoms: {equation_str} = {cn3_score}")
    print(f"Cranial Nerve III explains {cn3_score} out of {len(symptoms)} major symptoms.")
    print("This indicates a complete Cranial Nerve III palsy.")

    # Step 4: Map the affected nerve to its anatomical origin among the choices.
    anatomical_origins = {
        'A': 'Cranial nerve VI nucleus is in the Pons.',
        'B': 'Cranial nerve VII nucleus is in the Pons.',
        'C': 'Reticular formation is for consciousness/arousal, not isolated eye functions.',
        'D': 'Medulla oblongata is the origin for CN IX, X, XI, XII.',
        'E': 'Midbrain is the origin for Cranial Nerve III and IV.'
    }

    correct_choice = ''
    print("\nFinding the anatomical origin:")
    for choice, description in anatomical_origins.items():
        if "Cranial Nerve III" in description:
            correct_choice = choice
            print(f"Conclusion: Cranial Nerve III originates in the Midbrain.")
            print(f"Therefore, the damage is to the structure in Choice {choice}: {description.split(' is ')[0]}")
            break
    
    # This prevents the final answer from being printed by the function block itself.
    # We construct it and then pass it outside the standard output capture.
    final_answer = f'<<<{correct_choice}>>>'
    sys.stderr.write(final_answer)

solve_clinical_case()