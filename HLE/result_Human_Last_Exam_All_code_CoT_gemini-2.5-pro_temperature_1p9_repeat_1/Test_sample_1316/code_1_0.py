def bridges_experiment_simulation():
    """
    Simulates Bridges' experiment on non-disjunction in Drosophila
    to explain the origin of an X0 male.
    """

    # For simplicity, we'll model the classic white-eyed female x red-eyed male cross.
    # The X0 male gets his red-eye trait from his father's X chromosome.
    # Female genotype (white-eyed)
    female_genotype = "XwXw"
    # Male genotype (red-eyed)
    male_genotype = "X(R)Y"

    print("--- Simulating Bridges' Non-disjunction Experiment ---")
    print(f"Parental Cross: White-eyed Female ({female_genotype}) x Red-eyed Male ({male_genotype})\n")

    print("Step 1: Male Gamete Formation (Normal Meiosis)")
    male_gametes = ["X(R)", "Y"]
    print(f"The male produces two types of sperm: {male_gametes[0]} and {male_gametes[1]}\n")

    print("Step 2: Female Gamete Formation (with Non-disjunction in Meiosis I)")
    print("In Meiosis I, the two homologous X chromosomes (Xw and Xw) fail to separate.")
    print("This leads to the formation of two types of abnormal eggs:")
    female_abnormal_gametes = ["XwXw", "0"]
    print(f"  - One egg gets both X chromosomes: '{female_abnormal_gametes[0]}'")
    print(f"  - One egg gets no sex chromosome: '{female_abnormal_gametes[1]}' (This is a nullo-X gamete)\n")

    print("Step 3: Fertilization with Abnormal Female Gametes")
    print("This results in four possible exceptional offspring:\n")

    # This loop demonstrates the Punnett square for the exceptional cross
    for f_gamete in female_abnormal_gametes:
        for m_gamete in male_gametes:
            offspring_genotype = f_gamete.replace('0', '') + m_gamete
            
            print(f"  Female Egg ({f_gamete}) + Male Sperm ({m_gamete}) --> Offspring Genotype: {offspring_genotype}")

            if offspring_genotype == "XwXwY":
                print("    - Phenotype: White-eyed Female (fertile)")
            elif offspring_genotype == "X(R)XwXw":
                print("    - Phenotype: Superfemale (usually does not survive)")
            elif offspring_genotype == "Y":
                 print("    - Phenotype: Y0 (lethal, does not survive)")
            elif offspring_genotype == "X(R)":
                # Add the '0' for clarity in the final explanation
                final_genotype = "X(R)0"
                print(f"    - Phenotype: Red-eyed Male (sterile)")
                print("    *** This is the unexpected male offspring described in the question! ***")
            print("-" * 20)

    print("\nConclusion:")
    print("The unexpected red-eyed male with an X0 chromosomal makeup arises from a specific fertilization event:")
    # Explicitly print the final equation
    print("The final equation is: Egg('0') + Sperm('X(R)') --> Offspring('X(R)0')")
    print("This event is only possible because of non-disjunction of the X chromosome during female meiosis.")

bridges_experiment_simulation()