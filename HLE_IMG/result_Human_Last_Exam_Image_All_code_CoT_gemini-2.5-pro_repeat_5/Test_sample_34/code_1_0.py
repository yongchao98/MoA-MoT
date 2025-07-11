def solve_kinship_problem():
    """
    Analyzes Lévi-Strauss kinship diagram and evaluates which societies fit the pattern.
    """
    # 1. Define the relationship structure from the diagram
    # + indicates a familiar/warm relationship
    # - indicates a formal/distant relationship
    diagram_signature = {
        "Father/Son": "-",
        "Mother's Brother/Sister's Son": "+",
        "Husband/Wife": "+",
        "Brother/Sister": "-",
    }

    # 2. Define the relationship structures for each society based on ethnographic data
    societies = {
        "Trobriand": {
            "type": "Matrilineal",
            "Father/Son": "+",
            "Mother's Brother/Sister's Son": "-",
            "Husband/Wife": "+",
            "Brother/Sister": "-",
        },
        "Siuoi": {
            "type": "Matrilineal",
            "Father/Son": "+",
            "Mother's Brother/Sister's Son": "-",
            "Husband/Wife": "+",
            "Brother/Sister": "-",
        },
        "Lake Kubutu": {
            "type": "Patrilineal",
            "Father/Son": "-",  # Typical for patrilineal systems (authority)
            "Mother's Brother/Sister's Son": "+",  # Typical for patrilineal systems (familiarity)
            "Husband/Wife": "Unknown",
            "Brother/Sister": "Unknown",
        },
        "Tonga": {
            "type": "Patrilineal",
            "Father/Son": "-",
            "Mother's Brother/Sister's Son": "+",
            "Husband/Wife": "+",
            "Brother/Sister": "-",
        },
        "Cherkess": {
            "type": "Patrilineal",
            "Father/Son": "-",
            "Mother's Brother/Sister's Son": "+",
            "Husband/Wife": "-",  # Contradicts the diagram
            "Brother/Sister": "-",
        }
    }

    print("Step 1: Analyze the Diagram")
    print("The diagram represents a system with the following relationship dynamics:")
    for rel, nature in diagram_signature.items():
        state = "Formal/Distant" if nature == "-" else "Familiar/Warm"
        print(f"- {rel}: {state} ({nature})")
    print("\nThis core pattern (Formal F/S, Familiar MB/ZS) is typical of certain patrilineal societies.\n")

    print("Step 2: Evaluate Each Society")
    analysis_results = {}
    for name, data in societies.items():
        # Check core dynamic match (F/S and MB/ZS)
        core_match = (data["Father/Son"] == diagram_signature["Father/Son"] and
                      data["Mother's Brother/Sister's Son"] == diagram_signature["Mother's Brother/Sister's Son"])
        
        # Check for a perfect match on all known points
        mismatches = [rel for rel, sign in diagram_signature.items()
                      if data[rel] != "Unknown" and data[rel] != sign]
        
        perfect_match = core_match and not mismatches

        if perfect_match:
            result = "Perfect Match"
        elif core_match:
            if mismatches:
                result = f"Fits core dynamic, but mismatches on {mismatches}"
            else:
                result = "Fits core dynamic (Patrilineal Pattern)"
        else:
            result = "Mismatch on core dynamic"

        analysis_results[name] = result
        print(f"- {name} ({data['type']}): {result}")

    print("\nStep 3: Evaluate Answer Choices")
    answer_choices = {
        "A": ["Trobriand", "Siuoi"],
        "B": ["Siuoi", "Lake Kubutu"],
        "C": ["Lake Kubutu", "Tonga"],
        "D": ["Tonga", "Cherkess"],
        "E": ["Cherkess", "Trobriand"],
    }
    
    print(f"A. {analysis_results['Trobriand']} & {analysis_results['Siuoi']}")
    print(f"B. {analysis_results['Siuoi']} & {analysis_results['Lake Kubutu']}")
    print(f"C. {analysis_results['Lake Kubutu']} & {analysis_results['Tonga']}")
    print(f"D. {analysis_results['Tonga']} & {analysis_results['Cherkess']}")
    print(f"E. {analysis_results['Cherkess']} & {analysis_results['Trobriand']}")

    print("\nConclusion: Choice C provides the best fit, as both Lake Kubutu and Tonga are patrilineal systems correctly represented by the diagram's core logic, and Tonga is a perfect match for all depicted relationships.")

solve_kinship_problem()