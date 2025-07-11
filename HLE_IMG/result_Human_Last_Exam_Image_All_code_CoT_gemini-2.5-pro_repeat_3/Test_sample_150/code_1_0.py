def solve_reaction_sequence():
    """
    This function tracks the transformations in the provided reaction scheme
    and identifies the name of each intermediate and the final product.
    """
    
    # Step 1: Friedel-Crafts Acylation
    # Benzene reacts with propanoyl chloride in the presence of AlCl3.
    # An acyl group is attached to the benzene ring.
    intermediate_1 = "1-phenylpropan-1-one"
    
    # Step 2: Electrophilic Aromatic Bromination
    # The acyl group on 1-phenylpropan-1-one is a meta-director.
    # Bromine is added to the meta position of the phenyl ring.
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    
    # Step 3: Catalytic Hydrogenation
    # H2/Pd reduces both the benzylic ketone to an alkane and the aryl halide to a C-H bond.
    intermediate_3 = "propylbenzene"
    
    # Step 4: Benzylic Bromination
    # NBS with a radical initiator (benzoyl peroxide) brominates the benzylic position of propylbenzene.
    # The benzylic carbon is the one directly attached to the phenyl ring.
    final_product = "(1-bromopropyl)benzene"
    
    # Print the name of each identified molecule
    print(f"Intermediate-1: {intermediate_1}")
    print(f"Intermediate-2: {intermediate_2}")
    print(f"Intermediate-3: {intermediate_3}")
    print(f"Final Product: {final_product}")

solve_reaction_sequence()