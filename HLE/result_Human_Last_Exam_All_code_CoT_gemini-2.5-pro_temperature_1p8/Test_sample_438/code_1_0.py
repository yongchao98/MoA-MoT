def assess_disease_risk(fold_increase):
    """
    Analyzes the effect of increased antigen presentation on disease risk.

    Args:
      fold_increase: The factor by which the presentation of a
                     disease-causing self-antigen is increased.
    """

    # We can represent the base level of presentation with a relative value of 1.
    base_presentation_level = 1

    # The problem states the new variant increases this by a factor of 1000.
    # The new level of presentation is calculated as follows:
    new_presentation_level = base_presentation_level * fold_increase

    print(f"Base relative presentation level: {base_presentation_level}")
    print(f"Fold increase due to HLA variant: {fold_increase}")
    print("\nCalculating the new presentation level:")
    print(f"New Presentation Level = {base_presentation_level} * {fold_increase} = {new_presentation_level}")

    print("\n--- Conclusion ---")
    if new_presentation_level > base_presentation_level:
        conclusion = "INCREASE"
        reason = (
            "A massive increase in the presentation of a disease-causing self-antigen "
            "dramatically raises the probability of activating autoreactive T-cells, "
            "which initiates the autoimmune disease."
        )
    elif new_presentation_level < base_presentation_level:
        conclusion = "DECREASE"
        reason = (
            "A decrease in the presentation of a disease-causing self-antigen "
            "would make it less likely for autoreactive T-cells to be activated."
        )
    else:
        conclusion = "NOT CHANGE"
        reason = "The level of antigen presentation remains the same."

    print(f"This would likely {conclusion} a person's risk of developing the disease.")
    print(f"Reason: {reason}")

# The fold increase given in the problem is 1000.
assess_disease_risk(1000)