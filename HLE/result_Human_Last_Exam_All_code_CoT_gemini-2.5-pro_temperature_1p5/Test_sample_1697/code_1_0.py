def get_reaction_product():
    """
    This function determines the product of a two-step organic reaction.
    Step 1: N,N-diethyl-3-dimethylaminobenzamide reacts with sec-BuLi/TMEDA.
            This is a directed ortho-metalation. The most acidic proton is at C-2,
            as it is ortho to both the strong directing diethylamide group and the
            dimethylamino group. Lithiation occurs at this position.
    Step 2: The resulting organolithium intermediate is quenched with methyl iodide (CH3I).
            The nucleophilic lithiated carbon attacks the electrophilic methyl group,
            substituting the lithium with a methyl group.
    """
    starting_material = "N,N-diethyl-3-dimethylaminobenzamide"
    # The reaction adds a methyl group at position 2.
    final_product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    print(f"The final compound obtained is: {final_product_name}")

get_reaction_product()