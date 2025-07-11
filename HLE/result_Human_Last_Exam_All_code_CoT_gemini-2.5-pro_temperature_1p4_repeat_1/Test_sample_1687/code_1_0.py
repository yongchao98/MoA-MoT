def diagnose_patient():
    """
    This script evaluates potential diagnoses based on clinical findings
    from a patient's case after a colonoscopy.
    It uses a scoring system where points are awarded based on how
    well a diagnosis explains the key symptoms and facts.
    """

    diagnoses = {
        "A. Colonic perforation": {
            "description": "A hole in the wall of the colon.",
            "points": [1, 1],
            "rationale": "Points for being a colonoscopy complication (1) and causing internal bleeding (1), but the pain location is less typical."
        },
        "B. Lower GI bleeding": {
            "description": "Bleeding into the lower gastrointestinal tract.",
            "points": [1],
            "rationale": "Points for bleeding post-colonoscopy (1), but it doesn't explain LUQ/shoulder pain or peritoneal signs."
        },
        "C. Splenic laceration": {
            "description": "A tear in the spleen, a highly vascular organ in the LUQ.",
            "points": [3, 3, 1],
            "rationale": "High points for perfectly explaining the LUQ + left shoulder pain (3), massive hemorrhage (3), and being a known complication of a difficult colonoscopy (1)."
        },
        "D. Postpolypectomy syndrome": {
            "description": "Inflammation after removal of a polyp.",
            "points": [-99],
            "rationale": "Disqualified because the case explicitly states 'No polypectomy was performed' (-99)."
        }
    }

    print("Evaluating diagnoses based on key clinical findings...\n")

    highest_score = -float('inf')
    most_likely_diagnosis_letter = ''
    most_likely_diagnosis_name = ''

    # The letter 'C' from the choices corresponds to index 2
    # The letter from the user prompt corresponds to the following index: A=0, B=1, C=2, D=3
    final_choice_letter = "C"

    for i, (name, data) in enumerate(diagnoses.items()):
        score = sum(data["points"])
        
        # Build the equation string
        equation_parts = [str(p) for p in data["points"]]
        equation_str = " + ".join(equation_parts)

        print(f"Diagnosis: {name}")
        print(f"   Rationale: {data['rationale']}")
        print(f"   Score Calculation: {equation_str} = {score}")
        print("-" * 20)

        if score > highest_score:
            highest_score = score
            # Extracts the letter, e.g., "A" from "A. Colonic perforation"
            most_likely_diagnosis_letter = name.split('.')[0]
            most_likely_diagnosis_name = name


    print(f"\nConclusion: The diagnosis with the highest score is '{most_likely_diagnosis_name}' with a score of {highest_score}.")
    print("This is because it is the only condition that accounts for the combination of a difficult colonoscopy, left upper quadrant pain, left shoulder pain (Kehr's sign), and massive intra-abdominal hemorrhage leading to shock.")

    # The final answer must be in the specified format
    print(f"\nFinal Answer: The most likely diagnosis corresponds to choice {most_likely_diagnosis_letter}.")
    # Do not remove the line below, it's the final answer format.
    print(f'<<<{final_choice_letter}>>>')


diagnose_patient()