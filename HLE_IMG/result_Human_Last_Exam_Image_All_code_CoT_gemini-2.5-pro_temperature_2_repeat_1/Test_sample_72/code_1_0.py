def explain_pericyclic_reactions():
    """
    This function explains the two pericyclic reactions occurring in the given transformation.
    """
    explanation = [
        "The overall transformation occurs via two sequential photochemically allowed pericyclic reactions.",
        "",
        "Step 1: Electrocyclic Reaction",
        "The first step is the photochemical valence isomerization of hexafluorobenzene to form hexafluorobicyclo[2.2.0]hexa-2,5-diene (hexafluoro Dewar benzene).",
        "This is an example of an electrocyclic reaction, where a conjugated π system undergoes a ring closure.",
        "Reaction: Hexafluorobenzene --(hv)--> Hexafluoro Dewar benzene",
        "",
        "Step 2: [2+2] Cycloaddition",
        "The second step is the reaction between the intermediate, hexafluoro Dewar benzene, and cyclobutene.",
        "One of the 2π-electron systems (a double bond) of the Dewar benzene reacts with the 2π-electron system of cyclobutene.",
        "This forms a new four-membered ring, characteristic of a [2+2] cycloaddition.",
        "Reaction: Hexafluoro Dewar benzene + Cyclobutene --(hv)--> Final Product",
        "",
        "Conclusion:",
        "The two types of pericyclic reactions involved are an Electrocyclic reaction and a [2+2] cycloaddition."
    ]
    
    print('\n'.join(explanation))

explain_pericyclic_reactions()