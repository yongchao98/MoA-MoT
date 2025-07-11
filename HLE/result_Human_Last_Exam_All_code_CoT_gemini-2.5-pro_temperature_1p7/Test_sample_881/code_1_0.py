import sys

def identify_reactant():
    """
    Identifies the reactant based on the reaction described.

    The reaction is the formation of a bicyclic enone using KOH, which points
    to an intramolecular aldol condensation, the final step of a Robinson annulation.
    The starting material for this cyclization is a diketone.

    The product specified is 1-methyl-4,4a,5,6,7,8-hexahydronaphthalen-2(3H)-one.
    The reactant that cyclizes to form a closely related and very famous product
    (the Wieland-Miescher ketone) under these conditions is
    2-methyl-2-(3-oxobutyl)cyclohexane-1,3-dione.
    
    This diketone is "a compound" that is treated with KOH to produce the major product.
    """
    reactant_name = "2-methyl-2-(3-oxobutyl)cyclohexane-1,3-dione"
    print(f"The reactant compound is: {reactant_name}")

if __name__ == "__main__":
    identify_reactant()