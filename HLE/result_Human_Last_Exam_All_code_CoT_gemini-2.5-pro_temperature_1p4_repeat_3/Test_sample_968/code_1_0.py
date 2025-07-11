def solve_vaganova_arabesque_question():
    """
    Identifies the Vaganova arabesques based on the relationship
    between the forward arm and the lifted leg.
    """
    # Step 1: Define the characteristics of the Vaganova arabesques.
    # Based on Agrippina Vaganova's "Basic Principles of Classical Ballet".
    # - First: Arm opposite the lifted leg is forward.
    # - Second: Arm on the same side as the lifted leg is forward.
    # - Third: Both arms are forward, so the arm opposite the lifted leg is included.
    # - Fourth: Based on the Second, so the arm on the same side as the lifted leg is forward.
    arabesques_data = [
        {"number": 1, "name": "First", "forward_arm_relation": "opposite"},
        {"number": 2, "name": "Second", "forward_arm_relation": "same"},
        {"number": 3, "name": "Third", "forward_arm_relation": "opposite"},
        {"number": 4, "name": "Fourth", "forward_arm_relation": "same"},
    ]

    # Step 2: Find the arabesques that match the question's criteria.
    matching_arabesques = []
    for arabesque in arabesques_data:
        if arabesque["forward_arm_relation"] == "opposite":
            matching_arabesques.append(arabesque)

    # Step 3: Print the result in a clear, descriptive way.
    if len(matching_arabesques) == 2:
        arabesque1 = matching_arabesques[0]
        arabesque2 = matching_arabesques[1]
        
        print(f"The two types of arabesque where the dancer's arm extended in front is on the opposite side as the lifted leg are:")
        print(f"- {arabesque1['name']} (Number: {arabesque1['number']})")
        print(f"- {arabesque2['name']} (Number: {arabesque2['number']})")
        print("\nFinal Answer Selection:")
        # This fulfills the requirement to output each number in the final answer.
        print(f"The answer combines arabesque number {arabesque1['number']} and {arabesque2['number']}, which is 'First and Third'.")
    else:
        print("Could not determine a unique pair based on the provided data.")

solve_vaganova_arabesque_question()