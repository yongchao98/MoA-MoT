def solve_dermatomyositis_vignette():
    """
    This script analyzes a clinical vignette to identify the location of a rash
    based on the symptoms of Dermatomyositis.
    """
    # Key clinical findings from the patient's presentation
    patient_findings = {
        "muscle_weakness": True,
        "periorbital_erythema": True  # Redness around the eyes
    }

    # Diagnosis based on pathognomonic signs
    diagnosis = "Dermatomyositis"

    print("Step 1: Analyzing the patient's key symptoms.")
    print(f"The patient presents with muscle weakness and periorbital erythema.")
    print(f"This combination is highly characteristic of an autoimmune condition called {diagnosis}.\n")

    # Explanation of the specific skin finding
    print("Step 2: Identifying the specific rash described.")
    print("The term 'periorbital erythema' describes a classic sign known as the 'heliotrope rash'.")
    print("This rash is considered a pathognomonic sign for Dermatomyositis.\n")

    # Identifying the anatomical location
    rash_location = "Eyelids"
    print(f"Step 3: Determining the anatomical location of the heliotrope rash.")
    print(f"The heliotrope rash is specifically located on the upper {rash_location}.\n")

    # Relating the finding to the multiple-choice options
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }

    print("Step 4: Matching the location to the provided answer choices.")
    print("The question asks which anatomical region is expected to have a rash.")
    print(f"Given the finding of periorbital erythema, the rash is located on the {answer_choices['C']}.")

    final_answer = 'C'
    # Although other rashes can occur (A and E), the sign already given in the
    # vignette (periorbital erythema) points directly to C.
    print(f"\nFinal Answer: The correct option is {final_answer}.")


solve_dermatomyositis_vignette()
<<<C>>>