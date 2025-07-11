def explain_reactions():
    """
    Explains the two pericyclic reactions involved in the thermal transformation.
    """
    explanation = (
        "The thermal transformation shown occurs via two sequential pericyclic reactions:\n\n"
        "1. A 4π conrotatory electrocyclic ring-opening: The four-membered cyclobutene "
        "ring of the starting material opens. This is a 4π-electron system, and under thermal "
        "conditions (Δ), the reaction is conrotatory.\n\n"
        "2. A 6π disrotatory electrocyclic ring-closure: The resulting ten-membered ring "
        "intermediate contains a 1,3,5-hexatriene system. This 6π-electron system undergoes "
        "a thermal ring-closure, which is disrotatory, to form the final bicyclic product."
    )
    print(explanation)

if __name__ == "__main__":
    explain_reactions()