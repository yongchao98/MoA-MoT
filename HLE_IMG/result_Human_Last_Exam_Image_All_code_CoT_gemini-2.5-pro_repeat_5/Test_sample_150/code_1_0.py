def identify_product():
    """
    This function tracks the transformations in the provided reaction scheme
    to identify the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    # Benzene + Propanoyl chloride -> Propiophenone
    reactant1 = "Benzene"
    reagent1 = "Propanoyl chloride / AlCl3"
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: {reactant1} reacts with {reagent1} to form Intermediate-1.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Bromination
    # Propiophenone is meta-directing.
    reactant2 = intermediate_1
    reagent2 = "Br2 / FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {reactant2} is brominated with {reagent2}.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Catalytic Hydrogenation
    # H2/Pd reduces the benzylic ketone to an alkane and removes the bromine.
    reactant3 = intermediate_2
    reagent3 = "H2 / Pd"
    intermediate_3 = "Propylbenzene"
    print(f"Step 3: {reactant3} is reduced with {reagent3}.")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Radical Benzylic Bromination
    # NBS brominates the benzylic position.
    reactant4 = intermediate_3
    reagent4 = "NBS / (PhCO2)2, CCl4"
    final_product = "(1-Bromopropyl)benzene"
    print(f"Step 4: {reactant4} undergoes benzylic bromination with {reagent4}.")
    print(f"The final product is: {final_product}\n")

if __name__ == "__main__":
    identify_product()