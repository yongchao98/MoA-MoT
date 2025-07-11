import textwrap

def find_correct_statement():
    """
    This function evaluates several statements about the plant protein HR4
    and prints the correct one with an explanation.
    """

    # Dictionary containing the statements, their correctness (True/False), and an explanation.
    options = {
        'A': ("It is an interactor of the actin assembly factor ADF3",
              False,
              "Incorrect. There is no established direct interaction between the TIR-NLR protein HR4 and the actin-related protein ADF3."),
        'B': ("It contributes to the defense against the broad spectrum of powdery mildew pathogens",
              False,
              "Incorrect. The opposite is true. Research shows that the loss of HR4 function enhances resistance, indicating HR4 is a negative regulator of this specific defense pathway."),
        'C': ("It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen",
              False,
              "Incorrect. While functionally associated with proteins that do target the EHM (like RPW8.2), the HR4 protein itself is not known to be localized to the extrahaustorial membrane."),
        'D': ("It regulates the defense modulator PAD4 in plant defense against the Psm.",
              True,
              "Correct. HR4 is a TIR-NLR immune receptor. In response to pathogens like Psm (Pseudomonas syringae pv. maculicola), TIR-NLRs become activated and function as enzymes to produce small-molecule signals. These signals are perceived by the EDS1-PAD4 complex, thus regulating its activity to mount a defense response."),
        'E': ("HR4 is one of the interactors of the PAD4",
              False,
              "Incorrect. This statement is less precise than D. The regulation of PAD4 by HR4 is indirect, mediated by a small molecule, not a direct protein-protein interaction.")
    }

    print("Analyzing which statement about HR4 is true:\n")

    correct_key = None
    for key, (statement, is_correct, explanation) in options.items():
        if is_correct:
            correct_key = key
            print(f"The correct statement is: [{key}] {statement}")
            print("\nExplanation:")
            # Use textwrap for better formatting of the explanation
            wrapper = textwrap.TextWrapper(width=80, initial_indent="  ", subsequent_indent="  ")
            print(wrapper.fill(explanation))
            break
            
    if correct_key is None:
        print("Could not determine the correct statement.")

if __name__ == '__main__':
    find_correct_statement()