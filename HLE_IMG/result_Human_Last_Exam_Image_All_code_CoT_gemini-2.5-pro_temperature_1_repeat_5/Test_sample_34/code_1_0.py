def solve_kinship_diagram():
    """
    Analyzes the Lévi-Strauss kinship diagram and determines which societies it represents.
    """
    # Step 1 & 2: Define the relationships shown in the diagram
    diagram_pattern = {
        "Brother-Sister": "+ (Positive)",
        "Husband-Wife": "- (Negative/Formal)",
        "Father-Son": "+ (Positive/Affectionate)",
        "Mother's Brother-Sister's Son": "- (Negative/Authoritative)"
    }

    # Step 3 & 4: Define the theoretical patterns for different descent systems
    # In matrilineal systems, authority flows from the mother's brother.
    matrilineal_theory = {
        "Father-Son": "+ (Positive/Affectionate)",
        "Mother's Brother-Sister's Son": "- (Negative/Authoritative)"
    }
    # In patrilineal systems, authority flows from the father.
    patrilineal_theory = {
        "Father-Son": "- (Negative/Authoritative)",
        "Mother's Brother-Sister's Son": "+ (Positive/Affectionate)"
    }

    # Step 5: Define the societies from the answer choices
    societies = {
        "Trobriand": {"system": "matrilineal", "pattern": matrilineal_theory},
        "Siuoi": {"system": "matrilineal", "pattern": matrilineal_theory},
        "Lake Kubutu": {"system": "patrilineal", "pattern": patrilineal_theory},
        "Tonga": {"system": "patrilineal", "pattern": patrilineal_theory},
        "Cherkess": {"system": "patrilineal", "pattern": patrilineal_theory}
    }

    # --- Print out the analysis ---
    print("### Analysis of the Kinship Diagram ###")
    print("\n1. Relationships interpreted from the diagram:")
    for relation, sign in diagram_pattern.items():
        print(f"- {relation}: {sign}")

    print("\n2. Comparing the diagram's core axis with anthropological models:")
    print(f"   - Diagram's Father-Son relationship is: {diagram_pattern['Father-Son']}")
    print(f"   - Diagram's Mother's Brother-Sister's Son relationship is: {diagram_pattern['Mother's Brother-Sister's Son']}")

    print("\n   This pattern matches the theoretical model for a MATRILINEAL system.")
    print("   The pattern is the inverse of the model for a PATRILINEAL system.")


    print("\n3. Evaluating the answer choices:")
    answer_choices = {
        "A": ["Trobriand", "Siuoi"],
        "B": ["Siuoi", "Lake Kubutu"],
        "C": ["Lake Kubutu", "Tonga"],
        "D": ["Tonga", "Cherkess"],
        "E": ["Cherkess", "Trobriand"]
    }

    correct_choice = None
    for choice, names in answer_choices.items():
        society1_name = names[0]
        society2_name = names[1]
        society1 = societies[society1_name]
        society2 = societies[society2_name]

        # A choice is correct if BOTH societies match the matrilineal pattern of the diagram
        is_match1 = society1['system'] == 'matrilineal'
        is_match2 = society2['system'] == 'matrilineal'
        
        print(f"\n- Choice {choice}: {society1_name} ({society1['system']}) and {society2_name} ({society2['system']})")
        if is_match1 and is_match2:
            print("  Result: CORRECT. Both are matrilineal systems that fit the diagram's structure.")
            correct_choice = choice
        else:
            print("  Result: INCORRECT. One or both systems do not fit the diagram.")

    print(f"\n### Conclusion ###")
    print(f"The diagram represents the structure found in matrilineal societies. Both the Trobriand and Siuoi are classic examples of such systems.")
    print(f"Therefore, the correct choice is {correct_choice}.")

solve_kinship_diagram()