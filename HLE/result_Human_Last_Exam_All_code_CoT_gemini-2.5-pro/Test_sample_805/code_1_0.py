def model_braveheart_expression():
    """
    This function models the expression level of the Braveheart gene
    during the differentiation of embryonic stem cells (ESCs) into heart cells.
    """
    print("Analyzing the expression pattern of the Braveheart gene...")
    print("-" * 60)

    # Step 1: In the undifferentiated Embryonic Stem Cell (ESC)
    # In this state, the cell is pluripotent and has not committed to a lineage.
    # Heart-specific genes are not expressed.
    esc_expression = 1  # Represents a very low, basal expression level
    print(f"Stage: Embryonic Stem Cell (undifferentiated)")
    print(f"Braveheart Expression Level: {esc_expression} (Low/Off)")
    print("-" * 60)

    # Step 2: As the ESC begins to differentiate into a heart cell
    # Braveheart is a master regulator for cardiac lineage commitment. Its expression
    # must be significantly upregulated to start the process. This represents the "increase".
    peak_expression = 100  # Represents a sharp increase and peak expression
    print("Stage: Early Differentiation (transition from ESC to heart cell)")
    print("To initiate heart development, Braveheart expression sharply increases.")
    # This corresponds to the first part of the correct answer choice.
    print(f"Expression Change Equation: {esc_expression} -> {peak_expression} (Increase)")
    print("-" * 60)

    # Step 3: In the later-stage differentiating heart cell
    # Once Braveheart has activated the downstream cardiac gene network, its own
    # expression is no longer needed at high levels and it decreases.
    late_diff_expression = 25 # Represents the decrease from the peak
    print("Stage: Later Differentiation (maturing heart cell)")
    print("After its peak, Braveheart expression decreases as the cell matures.")
    # This corresponds to the second part of the correct answer choice.
    print(f"Expression Change Equation: {peak_expression} -> {late_diff_expression} (Decrease)")
    print("-" * 60)

    print("Conclusion:")
    print("The expression pattern is a transient pulse: a sharp INCREASE as differentiation begins,")
    print("followed by a DECREASE as differentiation proceeds.")
    print("This matches choice C, interpreting 'in embryonic stem cells' as the transition out of the ESC state.")

model_braveheart_expression()