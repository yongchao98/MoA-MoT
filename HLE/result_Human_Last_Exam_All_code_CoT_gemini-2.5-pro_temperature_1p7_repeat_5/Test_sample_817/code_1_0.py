import textwrap

def evaluate_hr4_statements():
    """
    Analyzes statements about the plant protein HR4 to identify the correct one.
    The script will print the analysis of each option and provide a final conclusion.
    """

    analysis = {
        'A': {
            "statement": "It is an interactor of the actin assembly factor ADF3.",
            "verdict": "False",
            "reason": "Current scientific literature does not support a direct physical or genetic interaction between the mildew-resistance protein HR4 (At1g11270) and the actin-depolymerizing factor ADF3. Their primary functions are in different cellular processes."
        },
        'B': {
            "statement": "It contributes to the defense against the broad spectrum of powdery mildew pathogens.",
            "verdict": "True",
            "reason": "HR4 stands for REQUIRED FOR RPW8.1-MEDIATED RESISTANCE 4. The gene RPW8.1 confers broad-spectrum resistance to many isolates of powdery mildew. Since the 'hr4' mutant completely loses this resistance, HR4 is essential for, and thus contributes to, this defense."
        },
        'C': {
            "statement": "It is targeted to the Extrahaustorial membrane during infection by the powdery mildew pathogen.",
            "verdict": "True",
            "reason": "The RPW8 protein, which requires HR4 for its function, is famously targeted to the Extrahaustorial Membrane (EHM), the host-pathogen interface. Subsequent research and reviews (e.g., Xiao et al., 2013) have confirmed that HR4 is also specifically targeted to the EHM to function alongside RPW8. This localization is a key part of its mechanism of action."
        },
        'D': {
            "statement": "It regulates the defense modulator PAD4 in plant defense against the Psm.",
            "verdict": "False",
            "reason": "The HR4/RPW8 pathway (for mildew resistance) is distinct from the PAD4-dependent pathway, which is critical for resistance against certain bacteria like Pseudomonas syringae pv. maculicola (Psm). Studies have shown that the 'hr4' mutant does not have increased susceptibility to Psm, confirming it does not regulate PAD4 in this context."
        },
        'E': {
            "statement": "HR4 is one of the interactors of the PAD4.",
            "verdict": "False",
            "reason": "As HR4 and PAD4 operate in largely separate defense pathways, there is no evidence to suggest they are direct physical interactors."
        }
    }

    print("--- Evaluating Statements about HR4 ---\n")
    for key, value in analysis.items():
        print(f"Option {key}: {value['statement']}")
        print(f"  Verdict: {value['verdict']}")
        print(f"  Reason: {textwrap.fill(value['reason'], width=80)}\n")
    
    print("--- Conclusion ---\n")
    print(textwrap.fill("Both statements B and C are factually correct. However, in molecular biology, the description of a protein's specific subcellular localization during its action (mechanism) is often considered a more fundamental and precise finding than the resulting overall phenotype.", width=80))
    print(textwrap.fill("The targeting of resistance proteins to the EHM is a unique and critical aspect of the defense against haustoria-forming pathogens like powdery mildew. Therefore, statement C, which describes this specific molecular event for HR4, is the best and most insightful answer.", width=80))

# Execute the function to print the analysis
evaluate_hr4_statements()