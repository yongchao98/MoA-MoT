def solve_kinship_diagram():
    """
    Analyzes the Lévi-Strauss kinship diagram and identifies the systems it represents.
    """

    print("Step 1: Deconstructing the Lévi-Strauss Kinship Diagram")
    print("---------------------------------------------------------")
    print("The diagram uses the following symbols:")
    print("  Δ : Male")
    print("  o : Female")
    print("  = : Marriage (Alliance)")
    print("  ― : Sibling bond (Consanguinity)")
    print("  \\ : Filiation (Parent-Child)")
    print("  + : Positive, familiar, or affectionate relationship")
    print("  - : Negative, reserved, formal, or antagonistic relationship")
    print("\n")

    print("Step 2: Determining the Relationship Pattern in the Diagram")
    print("---------------------------------------------------------")
    print("Based on the symbols and signs, the four key relationships in the 'atom of kinship' are:")
    # Husband (Δ) = Wife (o) has a '-' sign.
    husband_wife = "-"
    # Sister (o) ― Brother (Δ) has a '+' sign.
    brother_sister = "+"
    # The filiation line (\) from the parents to the son (Δ) has a '+' sign.
    father_son = "+"
    # Lévi-Strauss's theory posits an inverse relationship between Father/Son and Mother's Brother/Sister's Son.
    # Since Father/Son is '+', Mother's Brother/Sister's Son must be '-'.
    mothers_brother_sisters_son = "-"
    
    print(f"  - Husband/Wife: Reserved ({husband_wife})")
    print(f"  - Brother/Sister: Familiar ({brother_sister})")
    print(f"  - Father/Son: Familiar ({father_son})")
    print(f"  - Mother's Brother/Sister's Son: Reserved ({mothers_brother_sisters_son})")
    print("\n")

    print("Step 3: Analyzing Canonical Kinship Systems")
    print("---------------------------------------------")
    print("Let's analyze the two main examples used by Lévi-Strauss:")
    
    # Trobriand System (Matrilineal)
    trobriand_fs = "+" # Father is a loving friend
    trobriand_mbss = "-" # Mother's Brother is the authority figure
    trobriand_bs = "-" # Brother/Sister relationship has a strong avoidance taboo
    trobriand_hw = "+" # Husband/Wife relationship is relatively familiar
    print("A. Trobriand (Matrilineal) System:")
    print(f"  - Father/Son: Familiar ({trobriand_fs})")
    print(f"  - Mother's Brother/Sister's Son: Reserved ({trobriand_mbss})")
    print(f"  - Brother/Sister: Reserved ({trobriand_bs})")
    print(f"  - Husband/Wife: Familiar ({trobriand_hw})")
    print("   Comparison: Matches the diagram's vertical axis (F/S vs MB/SS) but contradicts the horizontal axis (B/S vs H/W).\n")

    # Cherkess System (Patrilineal)
    cherkess_fs = "-" # Father is the stern authority figure
    cherkess_mbss = "+" # Mother's Brother is an indulgent, familiar figure
    cherkess_bs = "+" # Brother/Sister relationship is familiar and free
    cherkess_hw = "-" # Husband/Wife relationship is marked by extreme reserve
    print("B. Cherkess (Patrilineal) System:")
    print(f"  - Father/Son: Reserved ({cherkess_fs})")
    print(f"  - Mother's Brother/Sister's Son: Familiar ({cherkess_mbss})")
    print(f"  - Brother/Sister: Familiar ({cherkess_bs})")
    print(f"  - Husband/Wife: Reserved ({cherkess_hw})")
    print("   Comparison: Matches the diagram's horizontal axis (B/S vs H/W) but contradicts the vertical axis (F/S vs MB/SS).\n")
    
    print("Step 4: Conclusion")
    print("------------------")
    print("The diagram represents an idealized theoretical structure. Neither of the primary ethnographic examples (Trobriand, Cherkess) perfectly fits the specific combination of signs in the diagram.")
    print("The diagram's vertical axis (Familiar Father / Reserved Uncle) is characteristic of the Trobriand (matrilineal) model.")
    print("The diagram's horizontal axis (Familiar Sister / Reserved Wife) is characteristic of the Cherkess (patrilineal) model.")
    print("Therefore, the diagram is best understood as an abstraction within a theory that is defined by the structural opposition between these two types of systems.")
    print("The question asks which *systems* are correctly represented. In the context of Lévi-Strauss's theory, the most relevant systems are the two poles he uses to build the entire framework: Cherkess and Trobriand.")
    print("Thus, the answer identifies the two foundational systems that the theory and its diagrams are meant to analyze and contrast.")

solve_kinship_diagram()