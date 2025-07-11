def braveheart_expression():
    """
    This function models and explains the expression of the Braveheart gene
    in embryonic stem cells vs. differentiating heart cells.
    """
    # 1. Define relative expression levels based on biological function.
    # In embryonic stem cells (ESCs), Braveheart is not needed, so its expression is very low.
    expression_esc = 1.0

    # During differentiation into heart cells, Braveheart is activated to guide the process.
    # Its expression increases significantly.
    expression_differentiating_heart = 15.0

    print("Analyzing the expression of the Braveheart gene:")
    print("-" * 45)
    print(f"In undifferentiated Embryonic Stem Cells, the relative expression level is low: {expression_esc}")
    print(f"In differentiating Heart Cells, the relative expression level increases significantly: {expression_differentiating_heart}")
    print("-" * 45)
    print("This shows that as stem cells commit to becoming heart tissue, the expression of the Braveheart gene is upregulated.")

    # Per the instructions, create a final equation showing the change.
    change_in_expression = expression_differentiating_heart - expression_esc
    print("\nRepresenting the change with an equation:")
    print(f"Final Expression ({expression_differentiating_heart}) - Initial Expression ({expression_esc}) = Change ({change_in_expression})")

braveheart_expression()