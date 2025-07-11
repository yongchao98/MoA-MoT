def explain_kinship_diagram():
    """
    Analyzes the Lévi-Strauss kinship diagram and identifies the systems it represents.
    """
    
    # Step 1 & 2: Define the structure shown in the diagram
    print("1. Analysis of the Diagram's Structure")
    print("=======================================")
    print("The diagram shows four individuals and their relationships:")
    print("  - Symbols: Male (Δ), Female (o)")
    print("  - Relationships: Marriage (=), Sibling (horizontal line), Descent (diagonal line)")
    print("  - Attitudes: Familiar/Warm (+), Reserved/Hostile (-)")
    print("\nThe specific relationships and attitudes in the diagram are:")
    
    brother_sister_rel = '+'
    husband_wife_rel = '-'
    maternal_uncle_nephew_rel = '+'
    # The father/son relationship is deduced from the structural principle.
    father_son_rel = '-'

    print(f"  - Brother/Sister: Familiar ({brother_sister_rel})")
    print(f"  - Husband/Wife: Reserved ({husband_wife_rel})")
    print(f"  - Mother's Brother/Sister's Son (Avuncular): Familiar ({maternal_uncle_nephew_rel})")
    
    # Step 3: Explain the underlying theory and deduce the final relationship
    print("\nLévi-Strauss's theory of the 'atom of kinship' posits a structural balance.")
    print("The rule is: (Uncle/Nephew rel) is to (Brother/Sister rel) as (Father/Son rel) is to (Husband/Wife rel).")
    print(f"Applying the rule to the diagram: ({maternal_uncle_nephew_rel}) is to ({brother_sister_rel}) as (Father/Son) is to ({husband_wife_rel}).")
    print("For the structure to be balanced, the Father/Son relationship must be opposite to the Uncle/Nephew one, resulting in:")
    print(f"  - Father/Son: Reserved ({father_son_rel})")

    # Step 4: Evaluate the primary case studies
    print("\n2. Evaluating Anthropological Case Studies")
    print("==========================================")
    print("We now compare this diagrammatic structure to real-world kinship systems.")
    
    print("\nCase A: The Cherkess (Patrilineal)")
    print("  - In this system, descent and authority pass from father to son.")
    print(f"  - To balance the system, the father-son bond is formal and reserved ({father_son_rel}).")
    print(f"  - The relationship with the mother's brother (the avuncular relative) is warm and familiar ({maternal_uncle_nephew_rel}).")
    print("  - CONCLUSION: The Cherkess system is a PERFECT MATCH for the specific pattern in the diagram.")

    print("\nCase B: The Trobriand Islanders (Matrilineal)")
    print("  - In this system, descent and identity pass through the mother's line.")
    print(f"  - The mother's brother is the figure of authority and discipline, making the relationship reserved ('-').")
    print(f"  - The biological father is a warm, loving friend, making the relationship familiar ('+').")
    print("  - This creates the INVERSE pattern to the one in the diagram.")
    print("  - CONCLUSION: The Trobriand system follows the same structural principle of opposition, but in its inverse form.")
    
    # Step 5: Conclude based on the theoretical nature of the diagram
    print("\n3. Final Conclusion")
    print("===================")
    print("The question asks which systems are correctly represented by the THEORETICAL diagram.")
    print("The diagram illustrates a fundamental principle, not just one specific society.")
    print("This principle is best demonstrated by Lévi-Strauss's two key, opposing examples:")
    print("  1. The Cherkess, who match the diagram's pattern directly.")
    print("  2. The Trobrianders, who match the inverse pattern, thereby proving the theory's structural validity across different descent types.")
    print("\nTherefore, the pair that is best represented by the overall theory is Cherkess-patrilineal and Trobriand-matrilineal.")

explain_kinship_diagram()