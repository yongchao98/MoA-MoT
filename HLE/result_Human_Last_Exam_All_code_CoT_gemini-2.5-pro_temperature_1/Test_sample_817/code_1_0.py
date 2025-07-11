import textwrap

def explain_hr4_function():
    """
    This function provides an explanation for the correct answer regarding the protein HR4.
    """
    explanation = {
        'A': "Incorrect. HR4 is a ROP-interactive CRIB motif-containing protein. Its primary and defining interaction is with ROP GTPases, not directly with the actin assembly factor ADF3.",
        'B': "Incorrect. While HR4 is involved in defense against powdery mildew, its role is highly specialized at the host-pathogen interface (the EHM). Describing its contribution as providing 'broad spectrum' resistance is not the most accurate description of its known specific function.",
        'C': "Correct. Scientific studies have demonstrated that during an infection by biotrophic fungi like powdery mildew, HR4 is specifically recruited to the plant cell membrane that surrounds the fungal feeding structure (haustorium). This membrane is known as the Extrahaustorial membrane (EHM). This localization is critical for its function in plant defense.",
        'D': "Incorrect. PAD4 is a key modulator in salicylic acid-dependent defense pathways, often studied in the context of bacterial pathogens like Psm (Pseudomonas syringae pv. maculicola). HR4 is not known to be a primary regulator of PAD4.",
        'E': "Incorrect. Similar to D, while signaling pathways in plant immunity are interconnected, HR4 is not documented as a direct interactor of PAD4. Its key interactors are ROPs."
    }

    print("Analysis of Answer Choices:")
    for option, reason in explanation.items():
        print(f"Choice {option}:")
        print(textwrap.fill(reason, width=80))
        print("-" * 20)

    print("\nConclusion: The most accurate statement based on established research is C.")

# Execute the function to print the explanation.
explain_hr4_function()

print("\n<<<C>>>")