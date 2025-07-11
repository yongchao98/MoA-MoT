def identify_reagents():
    """
    Identifies and prints the chemical reagents for the given reaction scheme.
    """
    # Reagent A is identified from the conversion of compound 1 to 2.
    # The replacement of a C-O-C bridge with a C-N(NH2)-C bridge points to hydrazine.
    reagent_A = "Hydrazine (N2H4)"

    # Reagent B is identified from the conversion of compound 2 to 3.
    # The incorporation of an N-propyl group into the structure points to propylamine.
    reagent_B = "Propylamine (CH3CH2CH2NH2)"

    print(f"Reagent A is: {reagent_A}")
    print(f"Reagent B is: {reagent_B}")

identify_reagents()