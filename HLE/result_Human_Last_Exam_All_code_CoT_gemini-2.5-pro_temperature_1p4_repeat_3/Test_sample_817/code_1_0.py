import textwrap

def explain_hr4_function():
    """
    This function evaluates several statements about the plant protein HR4
    and prints an analysis to determine which one is true.
    """
    # Introduction to the analysis
    print("Analyzing the function and interactions of the plant protein HR4...\n")

    # Analysis of each option
    analysis = {
        'A': "It is an interactor of the actin assembly factor ADF3. There is no prominent evidence in the literature to support a direct interaction between HR4 and ADF3.",
        'B': "It contributes to the defense against the broad spectrum of powdery mildew pathogens. Research, such as that by Hao et al. (2016) in 'Frontiers in Plant Science', shows that knockout of the HR4 gene compromises Arabidopsis resistance to both host (G. cichoracearum) and non-host (B. graminis f. sp. hordei) powdery mildew fungi. This supports the claim of its role in defense against a broad spectrum of these pathogens.",
        'C': "It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen. The same study by Hao et al. (2016) demonstrated through fluorescence microscopy that HR4 protein is recruited to and accumulates at the Extrahaustorial membrane (EHM), which surrounds the fungal feeding structure (haustorium) inside the plant cell. This is a key mechanistic finding.",
        'D': "It regulates the defense modulator PAD4 in plant defense against the Psm. There is no direct evidence to suggest that HR4 regulates PAD4. PAD4 is a well-studied component of TIR-NB-LRR-mediated immunity, a pathway not typically associated with HR4.",
        'E': "HR4 is one of the interactors of the PAD4. There is no evidence from major protein interaction studies that HR4 directly interacts with PAD4."
    }

    # Print the analysis for each option
    for option, text in analysis.items():
        print(f"Option {option}:")
        print(textwrap.fill(text, width=80))
        print("-" * 20)

    # Conclusion
    conclusion = """
Conclusion:
Both statements B and C are strongly supported by the same key research paper. However, statement C describes a specific molecular mechanism—the localization of the protein to the direct site of pathogen interaction (the EHM). This localization is the reason *why* HR4 can contribute to defense (statement B). In molecular biology, a specific, demonstrated mechanism is often considered a more precise and fundamental truth. Therefore, statement C is the most accurate and specific description of HR4's role during powdery mildew infection.
"""
    print(conclusion)

explain_hr4_function()
print("<<<C>>>")