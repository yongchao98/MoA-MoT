import pandas as pd

def analyze_kinship():
    """
    Analyzes the Lévi-Strauss kinship diagram and evaluates which societal systems
    from the options are correctly represented by it.
    """
    # 1. Analyze the Diagram
    # The diagram shows a husband (father), wife (mother), wife's brother (maternal uncle), and their son.
    # '+' denotes a familiar/tender relationship.
    # '-' denotes a formal/antagonistic/authoritarian relationship.
    
    diagram_analysis = {
        "relationship": ["Father-Son", "Mother's Brother-Sister's Son", "Brother-Sister", "Husband-Wife (Implied)"],
        "diagram_sign": ["+", "-", "+", "-"]
    }
    
    print("Step 1: Analyzing the theoretical diagram")
    print("The diagram illustrates a specific configuration of relationships in the 'atom of kinship':")
    # Using an f-string for clear output.
    for i in range(len(diagram_analysis["relationship"])):
        rel = diagram_analysis["relationship"][i]
        sign = diagram_analysis["diagram_sign"][i]
        explanation = "familiar/tender" if sign == "+" else "formal/authoritarian"
        print(f"- {rel}: {sign} (indicating a {explanation} relationship)")
    print("\nThe Husband-Wife sign is inferred as '-' based on Lévi-Strauss's structural rule: (MB/ZS) * (H/W) = (F/S) * (B/Z)")
    print("Which becomes: (-) * (H/W) = (+) * (+), therefore H/W must be (-).\n")


    # 2. Evaluate Ethnographic Data for Each System
    systems_data = {
        "System": ["Trobriand (Matrilineal)", "Siuoi (Matrilineal)", "Lake Kubutu (Patrilineal)", "Tonga (Patrilineal)", "Cherkess (Patrilineal)"],
        "Father-Son": ["+", "+", "-", "-", "-"],
        "MB-ZS": ["-", "-", "+", "+", "+"],
        "Brother-Sister": ["-", "+", "N/A", "-", "+"],
        "Husband-Wife": ["+", "-", "N/A", "+", "-"]
    }
    
    systems_df = pd.DataFrame(systems_data)

    print("Step 2: Evaluating the real-world kinship systems")
    print("We compare the diagram's pattern to known ethnographic data for each system:")
    print(systems_df.to_string(index=False))
    print("\n")

    # 3. Compare and Conclude
    print("Step 3: Comparing the diagram to the systems")
    print("The primary principle in the diagram is the opposition between a familiar Father-Son (+) relationship and an authoritarian Mother's Brother-Sister's Son (-) relationship.")
    print("This pattern is characteristic of certain matrilineal systems.\n")
    print("Let's see which systems fit this core F/S (+) vs. MB/ZS (-) principle:")

    core_matches = []
    perfect_matches = []
    
    for index, row in systems_df.iterrows():
        system_name = row['System']
        # Check core principle
        if row['Father-Son'] == "+" and row['MB-ZS'] == "-":
            core_matches.append(system_name)
            # Check for perfect match across all signs
            if row['Brother-Sister'] == "+" and row['Husband-Wife'] == "-":
                perfect_matches.append(system_name)

    print(f"- Systems matching the core principle (F/S: +, MB/ZS: -): {', '.join(core_matches)}")
    print(f"- Systems perfectly matching all four relationships: {', '.join(perfect_matches)}\n")

    print("Conclusion:")
    print("The 'Siuoi-matrilineal' system is a perfect match for the diagram across all four relationships.")
    print("The 'Trobriand-matrilineal' system matches the core principle (the F/S vs MB/ZS opposition) which defines the structural type, even though its Brother-Sister relationship (-) differs from the diagram's specific representation (+).")
    print("The other systems (Lake Kubutu, Tonga, Cherkess) are patrilineal and show the inverse core principle (F/S: -, MB/ZS: +), so they do not fit.")
    print("\nGiven the options, the best answer includes the two systems that exemplify the diagram's core structural principle.")
    print("Therefore, the systems that can be correctly represented by the diagram, as a type, are Trobriand-matrilineal and Siuoi-matrilineal.")

analyze_kinship()