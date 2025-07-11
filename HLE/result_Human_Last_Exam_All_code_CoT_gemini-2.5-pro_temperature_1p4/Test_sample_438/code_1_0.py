import sys

def calculate_risk_change(fold_increase):
    """
    This function illustrates the impact of increased antigen presentation on disease risk.
    Although disease risk is a complex multifactorial trait, a dramatic increase
    in the presentation of a disease-causing self-antigen is a major contributing factor.

    Args:
        fold_increase (int): The factor by which antigen presentation is increased.
    """

    # We can think of the baseline probability of activating an autoreactive T-cell
    # as a very small number. Let's represent it as a base unit of 1 for simplicity.
    base_presentation_level = 1

    # The HLA variant increases this by a large factor.
    new_presentation_level = base_presentation_level * fold_increase

    print("This script conceptually models the change in antigen presentation.")
    print("-" * 60)
    print("Human Leukocyte Antigen (HLA) molecules present peptides to T-cells.")
    print("If an HLA variant presents a disease-causing self-antigen more efficiently,")
    print("it increases the chances of activating autoreactive T-cells, which can lead to autoimmunity.")
    print("\nLet's model the change:")
    print(f"Baseline presentation level (arbitrary units): {base_presentation_level}")
    print(f"Fold increase due to HLA variant: {fold_increase}")

    # The final code needs to output each number in the final equation.
    print("\nResulting Equation:")
    # sys.stdout.write is used to print without a newline for formatting the equation
    sys.stdout.write(f"{base_presentation_level} * {fold_increase} = {new_presentation_level}\n")

    print("\nConclusion:")
    print(f"The presentation of the disease-causing self-antigen is increased by {fold_increase}-fold.")
    print("This would dramatically increase the probability of initiating an autoimmune response.")
    print("Therefore, the person's risk of developing the disease would likely INCREASE.")
    print("-" * 60)


if __name__ == '__main__':
    # The problem states a 1000-fold increase.
    increase_factor = 1000
    calculate_risk_change(increase_factor)
