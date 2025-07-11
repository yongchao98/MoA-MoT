def identify_pericyclic_reactions():
    """
    Identifies and prints the two types of pericyclic reactions involved
    in the photochemical reaction of hexafluorobenzene and cyclobutene.
    """

    # The transformation is a known two-step photochemical process.
    
    # The first step is a valence isomerization of the aromatic ring.
    # This is a type of electrocyclic reaction.
    # In this case, it is a 6-pi electron system undergoing ring closure.
    reaction_1 = {
        "type": "Electrocyclic reaction",
        "pi_electrons": 6
    }

    # The second step is a cycloaddition between the non-aromatic intermediate
    # and the cyclobutene molecule. Each reactant contributes a 2-pi electron system.
    reaction_2 = {
        "type": "Cycloaddition",
        "pi_system": "[2+2]"
    }

    print("The transformation involves two sequential photochemically allowed pericyclic reactions:")
    print("-" * 70)
    
    # Using a f-string to output the numbers from the analysis.
    print(f"1. A {reaction_1['pi_electrons']}π {reaction_1['type']}: Hexafluorobenzene isomerizes to hexafluoro-Dewar-benzene.")
    
    # Using a f-string to output the numbers/characters from the analysis.
    print(f"2. A {reaction_2['pi_system']} {reaction_2['type']}: The hexafluoro-Dewar-benzene intermediate reacts with cyclobutene.")
    print("-" * 70)

if __name__ == '__main__':
    identify_pericyclic_reactions()