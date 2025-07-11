import textwrap

def explain_hr4_function():
    """
    Explains the function of the HR4 protein and identifies the correct statement.
    """
    explanation = """
    The protein HR4 (also known as NRIP2) is a key component in the plant immune system. To determine the correct statement, we analyze the options based on scientific literature:

    A. It is an interactor of the actin assembly factor ADF3: There is no strong evidence for a primary interaction between HR4 and ADF3.

    B. It contributes to the defense against the broad spectrum of powdery mildew pathogens: While HR4 is involved in defense against powdery mildew, option C is a more precise and accurate description of its mechanism.

    C. It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen: This statement is correct. Studies have demonstrated that during an infection by a powdery mildew fungus, the HR4 protein is specifically recruited to the Extrahaustorial Membrane (EHM). The EHM is the crucial interface between the plant cell and the fungus's feeding structure (the haustorium). This localization is critical for its role in the plant's defense response.

    D. It regulates the defense modulator PAD4 in plant defense against the Psm.: HR4 and PAD4 are both immune-related proteins, but they operate in different signaling pathways. HR4's primary role is not the regulation of PAD4.

    E. HR4 is one of the interactors of the PAD4: The primary interactors of HR4 are NLR-type immune receptors, not PAD4.
    """

    print("Step-by-step analysis of the options for HR4:")
    print(textwrap.dedent(explanation).strip())
    print("\n" + "="*50)
    print("Therefore, the correct statement is:")
    print("C. It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen")
    print("="*50)

if __name__ == "__main__":
    explain_hr4_function()