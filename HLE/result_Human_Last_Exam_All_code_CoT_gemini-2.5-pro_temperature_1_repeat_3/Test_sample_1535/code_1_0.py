def solve_medical_case():
    """
    This script analyzes clinical findings to determine the location of an expected rash.
    """
    # Step 1: Define the key clinical findings from the vignette
    patient_findings = {
        "Muscle Weakness": True,
        "Myalgia": True,  # Muscle pain
        "Periorbital Erythema": True,  # Redness around the eyes
    }

    print("Analyzing the patient's most distinctive signs and symptoms:")
    for finding, present in patient_findings.items():
        if present:
            print(f"- {finding}")

    # Step 2: Formulate a likely diagnosis based on the constellation of symptoms
    # The combination of muscle weakness and periorbital erythema is highly suggestive of Dermatomyositis.
    diagnosis = "Dermatomyositis"
    print(f"\nThe combination of muscle weakness and periorbital erythema strongly suggests a diagnosis of: {diagnosis}")

    # Step 3: Define the characteristic rashes of the suspected diagnosis
    dermatomyositis_rashes = {
        "Heliotrope Rash": "Eyelids",
        "Gottron's Papules": "Dorsum of the hands",
        "Shawl Sign": "Shoulders",
    }
    print(f"\n{diagnosis} has several characteristic rashes, including:")
    for rash, location in dermatomyositis_rashes.items():
        print(f"- {rash} (located on the {location})")

    # Step 4: Use the specific vignette clue to pinpoint the most likely answer.
    # The "equation" shows how the key finding leads to the answer.
    vignette_clue = "Periorbital Erythema"
    corresponding_rash = "Heliotrope Rash"
    rash_location = "Eyelids"
    answer_choice = "C"

    print("\nFinal Reasoning Process:")
    print(f"The key finding is '{vignette_clue}'.")
    print("Equation: " + f"'{vignette_clue}'" + " => " + f"'{corresponding_rash}'" + " => " + f"Location: '{rash_location}'")
    print(f"\nThe patient's 'periorbital erythema' is the clinical description for the {corresponding_rash}, which characteristically appears on the {rash_location}.")
    print(f"\nTherefore, the correct answer choice is {answer_choice}.")

solve_medical_case()