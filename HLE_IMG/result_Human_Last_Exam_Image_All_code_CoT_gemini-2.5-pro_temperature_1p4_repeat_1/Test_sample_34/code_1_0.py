import sys

def solve_kinship_diagram():
    """
    Analyzes a Lévi-Strauss kinship diagram and identifies the corresponding social systems.
    """
    # Step 1 & 2: Deconstruct the diagram and determine its kinship structure.
    # + indicates a familiar relationship; - indicates a formal/antagonistic one.
    # Diagram Analysis:
    # 1. Father(Δ) to Son(Δ): The descending line is marked with '+'. -> Father/Son = + (Familiar)
    # 2. Husband(Δ) and Wife(O): The marriage bond '=' is associated with a '-' over the husband. -> Husband/Wife = - (Formal)
    # 3. Sister(O) and Brother(Δ): The sibling line '―' has a '+' over the sister. -> Brother/Sister = + (Familiar)
    # 4. Mother's Brother(Δ) and Sister's Son(Δ): The MB has a '-' above him. In this system, this denotes his relationship
    #    with his nephew, which is inverse to the Fa/So relationship. -> Mother's Brother/Sister's Son = - (Formal)
    diagram_system = {
        "Father/Son": "+",
        "Mother's Brother/Sister's Son": "-",
        "Brother/Sister": "+",
        "Husband/Wife": "-"
    }

    # Step 3: Model the ethnographic data for each society.
    societies = {
        "Trobriand-matrilineal": {
            "notes": "Classic matrilineal case. Authority lies with the mother's brother, making the father an affectionate figure.",
            "Father/Son": "+", "Mother's Brother/Sister's Son": "-", "Brother/Sister": "+", "Husband/Wife": "-"
        },
        "Siuoi-matrilineal": {
            "notes": "A matrilineal system structurally similar to the Trobriand Islanders.",
            "Father/Son": "+", "Mother's Brother/Sister's Son": "-", "Brother/Sister": "+", "Husband/Wife": "-"
        },
        "Lake Kubutu-patrilineal": {
            "notes": "A typical patrilineal system where the father is the authority figure.",
            "Father/Son": "-", "Mother's Brother/Sister's Son": "+"
        },
        "Tonga-patrilineal": {
            "notes": "Patrilineal system where the father is a figure of respect and authority, and the maternal uncle is indulgent.",
            "Father/Son": "-", "Mother's Brother/Sister's Son": "+"
        },
        "Cherkess-patrilineal": {
            "notes": "Patrilineal system with extreme formality between father and son, and familiarity with the maternal uncle.",
            "Father/Son": "-", "Mother's Brother/Sister's Son": "+"
        }
    }

    # Answer choices
    choices = {
        "A": ["Trobriand-matrilineal", "Siuoi-matrilineal"],
        "B": ["Siuoi-matrilineal", "Lake Kubutu-patrilineal"],
        "C": ["Lake Kubutu-patrilineal", "Tonga-patrilineal"],
        "D": ["Tonga-patrilineal", "Cherkess-patrilineal"],
        "E": ["Cherkess-patrilineal", "Trobriand-matrilineal"]
    }

    # Function to check for a match
    def check_match(society_name, society_data, diagram_data):
        is_match = True
        print(f"  - Comparing {society_name}:")
        for relation, value in diagram_data.items():
            if relation in society_data and society_data[relation] != value:
                is_match = False
                print(f"    - Mismatch in {relation}: Diagram is '{value}', but society is '{society_data[relation]}'")
            elif relation in society_data:
                 print(f"    - Match in {relation}: Both are '{value}'")
        if is_match:
            print(f"  --> Result: {society_name} MATCHES the diagram's structure.")
        else:
            print(f"  --> Result: {society_name} DOES NOT MATCH the diagram's structure.")
        return is_match

    # Step 4 & 5: Compare, identify, and print the solution
    print("Step 1: Analyzing the Kinship Diagram")
    print("The diagram represents the following relationship pattern:")
    for relation, attitude in diagram_system.items():
        print(f"  - {relation}: {'Familiar' if attitude == '+' else 'Formal/Antagonistic'} ({attitude})")
    print("\n" + "="*50 + "\n")

    print("Step 2: Evaluating Each Answer Choice")
    correct_choice_letter = None
    for choice, society_pair in choices.items():
        print(f"\n--- Checking Choice {choice}: {society_pair[0]} and {society_pair[1]} ---")
        match1 = check_match(society_pair[0], societies[society_pair[0]], diagram_system)
        match2 = check_match(society_pair[1], societies[society_pair[1]], diagram_system)
        if match1 and match2:
            correct_choice_letter = choice
            print(f"\nCONCLUSION: Both societies in Choice {choice} match the diagram.")
            break
        else:
            print(f"\nCONCLUSION: Choice {choice} is incorrect.")

    if correct_choice_letter:
        print("\n" + "="*50 + "\n")
        print(f"The analysis confirms that the correct answer is Choice {correct_choice_letter}.")
        sys.stdout.write(f"<<<{correct_choice_letter}>>>")
    else:
        print("No correct choice found based on the provided data.")

solve_kinship_diagram()