import sys

def identify_pericyclic_reactions():
    """
    Identifies and describes the two pericyclic reactions in the given
    thermal transformation.
    """
    
    # Define the characteristics of the two sequential reactions
    reactions = [
        {
            "step": 1,
            "type": "Electrocyclic Ring-Opening",
            "pi_electrons": 4,
            "stereochemistry": "Conrotatory",
            "description": "The starting material undergoes a thermal electrocyclic opening of the 4-membered cyclobutene ring."
        },
        {
            "step": 2,
            "type": "Electrocyclic Ring-Closure",
            "pi_electrons": 6,
            "stereochemistry": "Disrotatory",
            "description": "The resulting 10-membered ring intermediate contains a 1,3,5-hexatriene system that undergoes thermal closure."
        }
    ]

    print("The thermal transformation proceeds through two sequential pericyclic reactions:\n")

    # Print details for each reaction
    for reaction in reactions:
        print(f"Reaction {reaction['step']}: {reaction['type']}")
        print("---------------------------------")
        print(f"Description: {reaction['description']}")
        print(f"Number of π electrons: {reaction['pi_electrons']}")
        print(f"Stereochemical Mode (Thermal): {reaction['stereochemistry']}")
        print("")

if __name__ == "__main__":
    identify_pericyclic_reactions()
