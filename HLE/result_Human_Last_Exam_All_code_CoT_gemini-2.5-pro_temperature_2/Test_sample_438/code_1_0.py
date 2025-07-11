import sys

def calculate_risk_change():
    """
    This function explains and calculates the effect of an HLA variant on antigen presentation.
    """
    # The initial probability of presenting a specific self-antigen is typically very low.
    # We'll use a hypothetical baseline value for demonstration.
    baseline_presentation_level = 1

    # The problem states the variant increases this probability by a specific factor.
    fold_increase = 1000

    # Calculate the new presentation level with the high-risk HLA variant.
    variant_presentation_level = baseline_presentation_level * fold_increase

    # Print the explanation and the result.
    print("This simulation explains the impact of an HLA variant on disease risk.")
    print("-" * 50)
    print("The Human Leukocyte Antigen (HLA) system presents peptides to immune cells.")
    print("For autoimmune diseases, the key event is the presentation of a 'self-antigen'.")
    print("\nLet's quantify the change:")
    print(f"Baseline Level of Self-Antigen Presentation: {baseline_presentation_level} (arbitrary units)")
    print(f"Increase Caused by the HLA Variant: {fold_increase}-fold")
    print("\nFinal Calculation:")
    # The prompt requests the numbers in the final equation be outputted.
    print(f"New Presentation Level = {baseline_presentation_level} * {fold_increase}")
    print(f"New Presentation Level with Variant = {variant_presentation_level} (arbitrary units)")
    print("-" * 50)

    print("\nConclusion:")
    print("This massive increase in the presentation of a disease-causing self-antigen means")
    print("that autoreactive immune cells are much more likely to be stimulated.")
    print("This stimulation leads to an immune attack on the body's own tissues.")
    print("\nTherefore, this variant would significantly INCREASE the person's risk of disease.")


# Run the function. Python 2/3 compatible print.
if __name__ == '__main__':
    try:
        # This function prints the result, there is no need to capture it.
        calculate_risk_change()
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
