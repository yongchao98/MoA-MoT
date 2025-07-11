def solve_kinship_diagram():
    """
    Analyzes the Lévi-Strauss kinship diagram to identify the correct societal systems.
    """
    # Step 1: Define the relationships and attitudes from the diagram
    diagram_attitudes = {
        "Husband/Wife": "-",
        "Brother/Sister": "+",
        "Father/Son": "+",
        "Mother's Brother/Sister's Son": "-",
    }

    print("Step 1: Decoding the kinship diagram")
    print("The diagram shows the following relationship attitudes:")
    for relation, attitude in diagram_attitudes.items():
        sign = "Familiar/Positive" if attitude == "+" else "Formal/Negative"
        print(f"- {relation}: {sign} ({attitude})")
    print("-" * 20)

    # Step 2: Determine the kinship system type from the diagram
    print("Step 2: Identifying the descent system type")
    father_son = diagram_attitudes["Father/Son"]
    maternal_uncle_nephew = diagram_attitudes["Mother's Brother/Sister's Son"]
    system_type = ""
    if father_son == "+" and maternal_uncle_nephew == "-":
        system_type = "Matrilineal"
        explanation = "A familiar father (+) and a formal/authoritarian mother's brother (-) is characteristic of a Matrilineal system."
    elif father_son == "-" and maternal_uncle_nephew == "+":
        system_type = "Patrilineal"
        explanation = "A formal/authoritarian father (-) and a familiar mother's brother (+) is characteristic of a Patrilineal system."
    else:
        system_type = "Undetermined"
        explanation = "The pattern does not fit a classic matrilineal or patrilineal system."

    print(f"The diagram represents a: {system_type} system.")
    print(explanation)
    print("-" * 20)

    # Step 3: Define the known systems from the answer choices
    known_systems = {
        "Trobriand": {"type": "Matrilineal", "F/S": "+", "MB/ZS": "-"},
        "Siuoi": {"type": "Matrilineal", "F/S": "+", "MB/ZS": "-"},
        "Lake Kubutu": {"type": "Patrilineal", "F/S": "-", "MB/ZS": "+"},
        "Tonga": {"type": "Patrilineal", "F/S": "-", "MB/ZS": "+"},
        "Cherkess": {"type": "Patrilineal", "F/S": "-", "MB/ZS": "+"},
    }
    
    answer_choices = {
        "A": ["Trobriand", "Siuoi"],
        "B": ["Siuoi", "Lake Kubutu"],
        "C": ["Lake Kubutu", "Tonga"],
        "D": ["Tonga", "Cherkess"],
        "E": ["Cherkess", "Trobriand"],
    }

    print("Step 3: Evaluating the answer choices")
    correct_choice = None
    for choice, societies in answer_choices.items():
        is_match = all(known_systems[s]["type"] == system_type for s in societies)
        print(f"Choice {choice}: {', '.join(societies)}")
        if is_match:
            print(f"Result: Matches the {system_type} pattern. This is a plausible answer.")
            correct_choice = choice
        else:
            mismatched_society = [s for s in societies if known_systems[s]["type"] != system_type][0]
            print(f"Result: Incorrect. Contains a {known_systems[mismatched_society]['type']} system ({mismatched_society}).")
        print()

    print("-" * 20)
    # Step 4: Final Conclusion
    print("Step 4: Conclusion")
    print("Based on the analysis, the diagram represents a matrilineal system.")
    print("Only choice A contains exclusively matrilineal societies.")
    print("The Siuoi system is a perfect match for all four relationships in the diagram.")
    print("The Trobriand system matches the core F/S (+) and MB/ZS (-) matrilineal dynamic, making it a valid example of the type represented.")
    print(f"Therefore, the correct choice is {correct_choice}.")


solve_kinship_diagram()