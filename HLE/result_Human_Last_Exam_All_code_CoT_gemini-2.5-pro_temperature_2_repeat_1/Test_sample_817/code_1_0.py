import textwrap

def analyze_hr4_statements():
    """
    Analyzes multiple-choice questions about the plant resistance protein HR4.
    """
    analysis = {
        'A': "False. There is no clear evidence in scientific literature for a direct physical interaction between the resistance protein HR4 and the actin assembly factor ADF3.",
        'B': "False. HR4 is known to confer resistance to downy mildew (*Hyaloperonospora arabidopsidis*), not the taxonomically distinct powdery mildew pathogens.",
        'C': "False. This statement is incorrect because HR4's function is not associated with powdery mildew. Its localization during infection by the wrong pathogen is not its established role.",
        'D': "True. HR4 is a TIR-NLR protein providing resistance against *Pseudomonas syringae pv. maculicola* (*Psm*). As an upstream sensor, its activation upon pathogen recognition initiates a signaling cascade that requires and activates the downstream modulator PAD4. In this context, HR4 is considered to regulate PAD4.",
        'E': "False. Current models of TIR-NLR signaling suggest that HR4 acts as an enzyme to produce a small-molecule signal that activates the EDS1-PAD4 complex. This is an indirect regulatory relationship, not a direct protein-protein interaction."
    }

    print("Analysis of statements regarding the plant protein HR4:\n")

    for option, explanation in analysis.items():
        # Use textwrap to format the explanation nicely
        wrapped_text = textwrap.fill(f"Choice {option}: {explanation}", width=80)
        print(wrapped_text)
        print("-" * 80)

    print("\nConclusion: Based on the analysis of published research, statement D is the most accurate description.")

if __name__ == "__main__":
    analyze_hr4_statements()