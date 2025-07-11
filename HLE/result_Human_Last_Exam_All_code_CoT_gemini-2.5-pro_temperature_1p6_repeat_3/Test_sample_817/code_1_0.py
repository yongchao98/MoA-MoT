import textwrap

def explain_answer():
    """
    Explains the reasoning behind choosing the correct answer about the HR4 protein.
    """
    explanation = """
    The question asks for a true statement about the plant protein HR4. Based on scientific literature:

    A. It is an interactor of the actin assembly factor ADF3: This is not a well-supported claim.

    B. It contributes to the defense against the broad spectrum of powdery mildew pathogens: This is true. The RPW8/HR gene family is known for providing broad-spectrum resistance to powdery mildew.

    C. It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen: This is also true. This specific localization to the membrane surrounding the fungal feeding structure (haustorium) is a key and defining feature of the HR4/RPW8 protein family's defense mechanism.

    D. It regulates the defense modulator PAD4 in plant defense against the Psm: This is an inaccurate description of their complex relationship.

    E. HR4 is one of the interactors of the PAD4: This is very likely true, as its close homolog RPW8.2 has been shown to interact directly with PAD4.

    Both B and C are well-established facts about HR4. However, statement C describes a specific and unique molecular mechanism—the *how* of its function at the site of infection—which is a hallmark characteristic of this protein family. Statement B describes the outcome of this mechanism. In molecular biology, a specific mechanistic detail like localization to a key interface is often considered a more precise answer. Therefore, C is the best choice.
    """
    print(textwrap.dedent(explanation).strip())
    print("\nCorrect Answer: C")

# Execute the function to provide the explanation.
explain_answer()