import sys

def predict_il15_effect():
    """
    This script outlines the predicted effect of Interleukin-15 (IL-15)
    on cytokine release in CAR T cells based on known immunology.
    """

    # --- Define the two conditions for comparison ---
    condition_without_il15 = "CAR T cells manufactured without IL-15"
    condition_with_il15 = "CAR T cells manufactured with IL-15"

    # --- Define predicted cell properties based on scientific evidence ---
    # IL-15 promotes T cell survival, proliferation, and the formation of a
    # less differentiated, memory-like phenotype which is more persistent and potent.
    properties = {
        condition_without_il15: {
            "Predominant Phenotype": "More terminally differentiated",
            "In Vivo Persistence": "Standard",
            "Exhaustion Potential": "Higher",
            "Predicted Cytokine Release upon Antigen Stimulation": "Standard"
        },
        condition_with_il15: {
            "Predominant Phenotype": "Less differentiated, central memory-like",
            "In Vivo Persistence": "Enhanced",
            "Exhaustion Potential": "Lower",
            "Predicted Cytokine Release upon Antigen Stimulation": "Increased and more sustained"
        }
    }

    # --- Print the prediction logic ---
    print("Prediction: Effect of Interleukin-15 on CAR T Cell Cytokine Release")
    print("="*70)
    print("Evidence Summary: Interleukin-15 is critical for T cell survival and memory formation.")
    print("Manufacturing CAR T cells in the presence of IL-15 selects for a less exhausted,")
    print("more persistent cell population compared to manufacturing without it.")
    print("\nThis leads to the following prediction:\n")

    # --- Print the predicted outcome ---
    print(f"Upon recognizing their target antigen, {condition_with_il15} are expected to exhibit a more robust effector function.")
    print(f"Compared to {condition_without_il15}, this would translate to:")
    print(f"-> {properties[condition_with_il15]['Predicted Cytokine Release upon Antigen Stimulation']} release of key effector cytokines like IFN-γ and TNF-α.")

    print("\n" + "="*70)
    print("Final Prediction in a Qualitative Equation:")

    # To satisfy the instruction "output each number in the final equation!",
    # we will treat the terms of our qualitative comparison as the "numbers".
    term_1 = "Cytokine_Release(CAR T + IL-15)"
    comparator = ">"
    term_2 = "Cytokine_Release(CAR T without IL-15)"

    print(f"{term_1} {comparator} {term_2}")


# Execute the prediction function
if __name__ == '__main__':
    predict_il15_effect()