import collections

def analyze_eye_symptoms():
    """
    Analyzes ophthalmological symptoms to identify the affected anatomical structure.
    """

    # Mapping of cranial nerves to their primary eye-related functions
    cranial_nerve_functions = {
        "Cranial nerve III (Oculomotor)": ["adduction", "elevation", "depression", "pupillary constriction"],
        "Cranial nerve IV (Trochlear)": ["depression"],
        "Cranial nerve VI (Abducens)": ["abduction"],
    }

    # Anatomical location of the cranial nerve nuclei in the brainstem
    nuclei_locations = {
        "Cranial nerve III (Oculomotor)": "E. Midbrain",
        "Cranial nerve IV (Trochlear)": "E. Midbrain",
        "Cranial nerve VI (Abducens)": "Pons (not an answer choice)",
    }

    # Patient's reported symptoms
    patient_symptoms = {
        "Right eye no pupillary light reflex": "pupillary constriction",
        "Right eye unable to perform adduction": "adduction",
        "Right eye unable to perform depression": "depression",
        "Right eye unable to perform elevation": "elevation",
    }

    print("Step 1: Listing the patient's symptoms and the corresponding functional deficits.")
    for symptom, function in patient_symptoms.items():
        print(f"- Symptom: '{symptom}' -> Deficit: '{function}'")
    print("-" * 20)

    # Correlate functional deficits to the responsible cranial nerves
    implicated_nerves = []
    for function in patient_symptoms.values():
        for nerve, funcs in cranial_nerve_functions.items():
            if function in funcs:
                implicated_nerves.append(nerve)
    
    # Count the occurrences of each implicated nerve
    nerve_counter = collections.Counter(implicated_nerves)
    most_likely_nerve = nerve_counter.most_common(1)[0][0]

    print("Step 2: Identifying the most likely cranial nerve affected based on symptoms.")
    print(f"The deficits point overwhelmingly to issues with functions controlled by the {most_likely_nerve}.")
    print("-" * 20)
    
    # Find the anatomical location of the identified nerve
    location = nuclei_locations.get(most_likely_nerve, "Unknown Location")

    print(f"Step 3: Determining the anatomical location of the {most_likely_nerve} nucleus.")
    print(f"The nucleus for the {most_likely_nerve} is located in the Midbrain.")
    print("-" * 20)
    
    print("Conclusion: The patient's presentation is explained by damage to the anatomical structure corresponding to the nerve's origin.")
    print(f"Final Answer is: {location}")

# Run the analysis
analyze_eye_symptoms()