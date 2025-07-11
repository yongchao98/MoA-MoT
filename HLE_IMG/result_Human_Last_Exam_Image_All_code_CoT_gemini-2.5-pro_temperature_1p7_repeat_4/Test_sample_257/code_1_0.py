import sys

def predict_nmr_signal():
    """
    This function outlines the reasoning for determining the splitting pattern and
    integration of the most deshielded proton in Compound 1.
    """
    
    # Analysis of the reaction and product
    reaction_type = "electrophilic aromatic sulfonation"
    product_feature = "two sulfonic acid (-SO3H) groups on the central aromatic ring"
    
    # Identification of the most deshielded proton
    proton_location = "the proton at position 5 on the central ring, between the two new -SO3H groups"
    
    # Determination of NMR properties
    splitting_pattern = "singlet"
    number_of_neighbors = 0
    integration = 1
    
    # Printing the explanation
    print("Chemical Analysis:")
    print(f"1. The reaction is an {reaction_type}.")
    print(f"2. The product, Compound 1, contains {product_feature}.")
    print(f"3. The most deshielded proton in Compound 1 is {proton_location}.")
    
    print("\nNMR Prediction:")
    print(f" - Splitting Pattern: The proton has {number_of_neighbors} adjacent proton neighbors.")
    print(f"   Therefore, according to the n+1 rule, its signal is a {splitting_pattern}.")
    
    # The prompt requests outputting each number in the final equation.
    # The integration value is the key number here.
    print(f" - Integration: There is only one such proton in the molecule.")
    print(f"   The integration value is: {integration}")
    
    print("\n---")
    print("Final Answer:")
    print(f"The highest deshielded proton peak has a splitting pattern of a {splitting_pattern} and an integration of {integration}H.")

# Execute the function to print the solution
predict_nmr_signal()
