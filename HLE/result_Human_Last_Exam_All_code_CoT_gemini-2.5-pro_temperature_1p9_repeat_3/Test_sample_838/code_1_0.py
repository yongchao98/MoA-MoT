def final_product_identity():
    """
    Prints the analysis of the final product C.
    """
    final_product = "The final product, C, is 1-bromo-4-phenylbutane."
    iupac_name = "IUPAC Name: 1-bromo-4-phenylbutane"
    chirality = ("Chirality Explanation: The final product is achiral. "
                 "The original stereocenter in [(3S)-3-bromobutyl]benzene was eliminated in the first reaction step to form 4-phenylbut-1-ene. "
                 "Since all subsequent intermediates were achiral, the final product, 1-bromo-4-phenylbutane, contains no chiral centers.")

    print("--- Final Product Analysis ---")
    print(final_product)
    print(iupac_name)
    print(chirality)

if __name__ == "__main__":
    final_product_identity()