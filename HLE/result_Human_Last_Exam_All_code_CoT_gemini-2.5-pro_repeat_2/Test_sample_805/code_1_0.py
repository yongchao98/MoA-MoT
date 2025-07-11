import sys

def analyze_braveheart_expression():
    """
    Analyzes and explains the expression pattern of the Braveheart (Bvht) lncRNA
    in embryonic stem cells and differentiating heart cells.
    """
    # Step 1: Define representative expression levels based on biological function.
    # In undifferentiated Embryonic Stem Cells (ESCs), lineage-specific regulators are off.
    esc_expression_level = 1  # Represents a very low, suppressed, or basal level.

    # During differentiation into heart cells, Braveheart is a key driver and is strongly activated.
    differentiating_heart_cell_expression = 95  # Represents a high, induced level of expression.

    # Step 2: Explain the biological context.
    print("### Analysis of Braveheart Gene Expression ###")
    print("\n1. Function of Braveheart (Bvht):")
    print("   Bvht is a crucial regulator that promotes the differentiation of embryonic stem cells (ESCs) into heart cells. It is essential for the formation of the cardiovascular lineage.")

    print("\n2. Expression in Different Cell States:")
    print(f"   - In undifferentiated ESCs, Bvht expression is kept very low (suppressed). Modeled level: {esc_expression_level}")
    print(f"   - To become a heart cell, Bvht expression must be significantly increased. Modeled level: {differentiating_heart_cell_expression}")

    # Step 3: Present the change as a simple equation with the modeled numbers.
    print("\n3. Calculating the Change:")
    print("   The change in expression when an ESC differentiates into a heart cell is an increase.")
    change_value = differentiating_heart_cell_expression - esc_expression_level
    print(f"   Equation: (Differentiating Cell Expression) - (ESC Expression) = Change")
    print(f"   Calculation: {differentiating_heart_cell_expression} - {esc_expression_level} = {change_value}")

    # Step 4: Evaluate the provided options.
    print("\n4. Conclusion:")
    print("   The expression of Braveheart is very low or suppressed in embryonic cells and it increases substantially in differentiating heart cells.")
    print("   Option D, 'It decreases in embryonic cells and slightly increases in differentiating heart cells,' is the best fit.")
    print("   - 'Decreases in embryonic cells' is likely imprecise wording for 'is low/suppressed'.")
    print("   - It correctly identifies the direction of change ('increases') during differentiation.")

# Execute the analysis
analyze_braveheart_expression()