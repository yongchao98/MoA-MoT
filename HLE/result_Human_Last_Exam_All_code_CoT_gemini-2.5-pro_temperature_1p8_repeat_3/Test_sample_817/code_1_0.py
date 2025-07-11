import textwrap

def analyze_hr4_statements():
    """
    Analyzes the multiple-choice options regarding the plant protein HR4 based on scientific literature.
    """
    # Introduction to HR4
    print("HR4 (HOP-LIKE-RESISTANCE4) is a Tudor-like Staphylococcal nuclease (SN) domain protein in Arabidopsis thaliana.")
    print("-" * 50)

    # Analysis of each option
    print("Analyzing the given statements:\n")

    # Option A
    analysis_A = "A. It is an interactor of the actin assembly factor ADF3. \n   - Verdict: False. There is no significant evidence linking HR4 (a Tudor-SN protein involved in nucleic acid interactions and immunity) with ADF3 (an actin-depolymerizing factor). They belong to different functional pathways."
    print(textwrap.fill(analysis_A, width=80))
    print("\n")

    # Option D
    analysis_D = "D. It regulates the defense modulator PAD4 in plant defense against the Psm. \n   - Verdict: False. Research shows that HR4's role is specific to powdery mildew defense. Mutants of HR4 do not show altered susceptibility to the bacterial pathogen Pseudomonas syringae (Psm), making a role in Psm defense unlikely."
    print(textwrap.fill(analysis_D, width=80))
    print("\n")

    # Option E
    analysis_E = "E. HR4 is one of the interactors of the PAD4. \n   - Verdict: False. There is no evidence for a direct interaction between HR4 and PAD4. HR4 is known to be required for the full function of another resistance protein, RPW8.2, in the context of powdery mildew defense."
    print(textwrap.fill(analysis_E, width=80))
    print("\n")
    
    # Options B and C are related and supported by the same key research.
    print("Now analyzing the most relevant options, B and C:\n")
    
    # Option B
    analysis_B = "B. It contributes to the defense against the broad spectrum of powdery mildew pathogens. \n   - Verdict: True. A key paper (Pennington et al., 2019, The Plant Cell) establishes that HR4 is essential for broad-spectrum resistance to powdery mildew pathogens."
    print(textwrap.fill(analysis_B, width=80))
    print("\n")

    # Option C
    analysis_C = "C. It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen. \n   - Verdict: True. The same paper shows that functional HR4 protein is specifically localized to the extrahaustorial membrane (EHM), the crucial interface between the plant cell and the invading fungal structure. This localization is key to its function."
    print(textwrap.fill(analysis_C, width=80))
    print("\n")

    # Conclusion
    print("-" * 50)
    conclusion = "Both B and C are factually correct. However, C provides a more specific, mechanistic detail about HR4's function. Its unique localization to the EHM explains *how* it contributes to defense. In molecular biology, this type of mechanistic insight is often considered a more precise and fundamental finding than the description of the overall phenotype. Therefore, C is the best answer."
    print(textwrap.fill(conclusion, width=80))


analyze_hr4_statements()