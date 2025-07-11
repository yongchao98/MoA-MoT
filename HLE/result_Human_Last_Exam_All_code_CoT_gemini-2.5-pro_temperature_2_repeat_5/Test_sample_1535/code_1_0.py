def find_rash_location():
    """
    This script analyzes the key clinical finding from the vignette to determine
    the most likely location of a rash.
    """

    # The most specific physical exam finding is "periorbital erythema".
    key_finding_term = "periorbital"
    finding_description = "erythema"

    # In medicine, "periorbital" refers to the area around the eyes.
    # A characteristic rash in this area is a key diagnostic sign.
    anatomical_location_map = {
        "periorbital": "Eyelids"
    }

    # Answer choices provided in the problem
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }

    # Determine the location based on the key term
    if key_finding_term in anatomical_location_map:
        expected_location = anatomical_location_map[key_finding_term]
        
        # Find the corresponding letter for the final answer
        correct_answer_letter = [letter for letter, location in answer_choices.items() if location == expected_location][0]

        print("The key physical finding is 'periorbital erythema'.")
        print("The term 'periorbital' points to the region around the eyes.")
        print("This suggests a heliotrope rash, which characteristically appears on the:")
        print(expected_location)
        print("This corresponds to answer choice:")
        print(correct_answer_letter)

find_rash_location()