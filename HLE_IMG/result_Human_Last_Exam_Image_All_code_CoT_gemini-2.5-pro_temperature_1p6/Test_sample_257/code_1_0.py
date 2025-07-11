import sys

def determine_nmr_properties():
    """
    Analyzes the sulfonation of Pr-DAOTA and determines the 1H NMR properties
    of the most deshielded proton in the resulting product (Compound 1).
    """

    # Step 1 & 2: Analyze the reaction and predict the product.
    # The reaction is an electrophilic aromatic sulfonation. The product, Compound 1,
    # is the disulfonated Pr-DAOTA cation, making it water-soluble.

    # Step 3: Identify the most deshielded proton.
    # The most deshielded proton in the triangulenium system is the one on the top ring,
    # positioned centrally between the two N-propylamino substituents. Its highly
    # electron-deficient environment leads to the highest chemical shift.

    # Step 4: Determine the integration and splitting pattern.

    # Integration: There is only one proton in this unique environment.
    integration = 1

    # Splitting Pattern: This proton is coupled to two equivalent adjacent (ortho) protons.
    # According to the n+1 rule, where n is the number of neighbors:
    # n = 2
    # Splitting = n + 1
    # This results in a triplet.
    splitting_pattern = "triplet"

    # Output the result.
    print("Based on the analysis of the structure of Compound 1:")
    print("The highest deshielded proton peak in the 1H NMR spectrum has:")
    print(f"- Splitting Pattern: {splitting_pattern}")
    # The instruction requires outputting the number from the "equation".
    # Here, the integration number is the key numerical result.
    print(f"- Integration: {integration}H")


# Execute the function to print the solution.
determine_nmr_properties()