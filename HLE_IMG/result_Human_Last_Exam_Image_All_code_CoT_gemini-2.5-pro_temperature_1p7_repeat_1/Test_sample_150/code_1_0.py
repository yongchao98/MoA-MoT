def solve_reaction_sequence():
    """
    This function tracks the transformations in the given chemical reaction sequence
    and identifies the final product.
    """
    
    # --- Step 1: Friedel-Crafts Acylation ---
    intermediate_1_name = "1-phenylpropan-1-one (Propiophenone)"
    
    # --- Step 2: Electrophilic Bromination ---
    intermediate_2_name = "1-(3-bromophenyl)propan-1-one"
    
    # --- Step 3: Catalytic Hydrogenolysis ---
    intermediate_3_name = "Propylbenzene"
    
    # --- Step 4: Radical Benzylic Bromination ---
    final_product_name = "1-bromo-1-phenylpropane"
    
    print("Tracking the reaction pathway:")
    print("-" * 35)
    
    print("Step 1: Benzene reacts with propanoyl chloride and AlCl3.")
    print(f"--> Intermediate-1: {intermediate_1_name}")
    print("-" * 35)
    
    print(f"Step 2: {intermediate_1_name} reacts with Br2/FeBr3.")
    print(f"--> Intermediate-2: {intermediate_2_name}")
    print("-" * 35)
    
    print(f"Step 3: {intermediate_2_name} reacts with H2/Pd.")
    print(f"--> Intermediate-3: {intermediate_3_name}")
    print("-" * 35)
    
    print(f"Step 4: {intermediate_3_name} reacts with NBS and a radical initiator.")
    print(f"--> Final Product: {final_product_name}")
    print("-" * 35)

# Execute the function to display the result.
solve_reaction_sequence()