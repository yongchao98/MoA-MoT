def get_reaction_product():
    """
    This function identifies the product of a directed ortho-metalation reaction.
    
    The reaction involves:
    1. Starting material: N,N-diethyl-3-dimethylaminobenzamide
    2. Reagents: 
        a) sec-BuLi/TMEDA for lithiation
        b) Methyl iodide (CH3I) as an electrophile
    
    The amide and amino groups direct the lithiation to the C2 position,
    which is then methylated by methyl iodide.
    """
    
    # Define parts of the molecule name
    amide_part = "N,N-diethyl"
    substituent_part = "2-methyl-3-dimethylamino"
    base_name = "benzamide"
    
    # Combine the parts to form the final product name
    final_product_name = f"{amide_part}-{substituent_part}{base_name}"
    
    print(f"The final product is: {final_product_name}")

get_reaction_product()