def solve_kinship_diagram():
    """
    Analyzes a Lévi-Strauss kinship diagram and identifies the corresponding societal systems.
    """

    # Step 1 & 2: Deconstruct the diagram and analyze the relationships
    print("Step 1: Analyzing the Lévi-Strauss Kinship Diagram")
    print("----------------------------------------------------")
    print("The diagram uses the following symbols:")
    print("  Δ: Male")
    print("  o: Female")
    print("  =: Marriage bond (Husband-Wife)")
    print("  ―: Sibling bond (Brother-Sister)")
    print("  Connecting lines from top to bottom imply descent (Parent-Child).")
    print("  +: A familiar, tender, or positive relationship.")
    print("  -: A formal, antagonistic, or negative relationship.")
    print("\nFrom the diagram, we can deduce the following four key relationships:")

    # Step 3: Establish the diagram's pattern
    diagram_pattern = {
        "Brother-Sister": "+",
        "Husband-Wife": "-",
        "Father-Son": "+",
        "Maternal Uncle-Nephew": "-"
    }

    print(f"  - Brother (right Δ) and Sister (o): The relationship is marked with a '{diagram_pattern['Brother-Sister']}' (Familiar).")
    print(f"  - Husband (left Δ) and Wife (o): The relationship is marked with a '{diagram_pattern['Husband-Wife']}' (Formal/Antagonistic).")
    print(f"  - Father (left Δ) and Son (bottom Δ): The relationship is marked with a '{diagram_pattern['Father-Son']}' (Familiar).")
    print(f"  - Mother's Brother (right Δ) and Sister's Son (bottom Δ): The relationship is marked with a '{diagram_pattern['Maternal Uncle-Nephew']}' (Formal/Antagonistic).")
    print("\nThis pattern is characteristic of systems where authority lies with the maternal uncle, not the father.")
    print("This is common in matrilineal societies.")
    print("----------------------------------------------------\n")

    # Step 4: Define the patterns for the ethnographic cases
    societies = {
        "Trobriand-matrilineal": {
            "description": "Classic matrilineal system. Authority resides with the mother's brother, who has a formal relationship with his nephew. The father has a warm, familiar relationship.",
            "pattern": {"Brother-Sister": "+", "Husband-Wife": "-", "Father-Son": "+", "Maternal Uncle-Nephew": "-"}
        },
        "Siuoi-matrilineal": {
            "description": "A matrilineal system similar to the Trobriand, where the maternal uncle holds authority.",
            "pattern": {"Brother-Sister": "+", "Husband-Wife": "-", "Father-Son": "+", "Maternal Uncle-Nephew": "-"}
        },
        "Lake Kubutu-patrilineal": {
            "description": "A patrilineal system where the father holds authority, leading to a formal father-son relationship.",
            "pattern": {"Father-Son": "-"} # Key differentiator
        },
        "Tonga-patrilineal": {
            "description": "A patrilineal system where the father's authority makes the father-son relationship formal.",
            "pattern": {"Father-Son": "-"} # Key differentiator
        },
        "Cherkess-patrilineal": {
            "description": "A patrilineal system that represents the inverse of the diagram. The father-son relationship is very formal, while the maternal uncle-nephew relationship is familiar and indulgent.",
            "pattern": {"Brother-Sister": "+", "Husband-Wife": "-", "Father-Son": "-", "Maternal Uncle-Nephew": "+"}
        }
    }

    # Step 5: Compare the diagram's pattern with each society
    print("Step 2: Comparing the Diagram's Pattern with Ethnographic Examples")
    print("----------------------------------------------------")

    def check_match(society_name, society_data):
        print(f"Checking: {society_name}")
        print(f"  - Description: {society_data['description']}")
        is_match = society_data.get("pattern") == diagram_pattern
        result = "MATCHES" if is_match else "DOES NOT MATCH"
        print(f"  - Result: {result} the diagram's pattern.\n")
        return is_match

    results = {name: check_match(name, data) for name, data in societies.items()}
    print("----------------------------------------------------\n")

    # Step 6: Evaluate the answer choices
    print("Step 3: Evaluating the Answer Choices")
    print("----------------------------------------------------")
    answer_choices = {
        "A": ["Trobriand-matrilineal", "Siuoi-matrilineal"],
        "B": ["Siuoi-matrilineal", "Lake Kubutu-patrilineal"],
        "C": ["Lake Kubutu-patrilineal", "Tonga-patrilineal"],
        "D": ["Tonga-patrilineal", "Cherkess-patrilineal"],
        "E": ["Cherkess-patrilineal", "Trobriand-matrilineal"]
    }

    correct_choice = None
    for choice, societies_in_choice in answer_choices.items():
        match1 = results.get(societies_in_choice[0], False)
        match2 = results.get(societies_in_choice[1], False)
        print(f"Choice {choice}: {societies_in_choice[0]} ({'Match' if match1 else 'No Match'}) and {societies_in_choice[1]} ({'Match' if match2 else 'No Match'})")
        if match1 and match2:
            correct_choice = choice
            print(f"Conclusion: Choice {choice} is correct because both systems fit the pattern shown in the diagram.")

    return correct_choice

final_answer = solve_kinship_diagram()
print(f"<<<{final_answer}>>>")
