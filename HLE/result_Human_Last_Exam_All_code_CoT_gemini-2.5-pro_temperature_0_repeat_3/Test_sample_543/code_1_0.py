import pubchempy as pcp

def find_iupac_name():
    """
    Finds the IUPAC name of the reaction product, o-terphenyl,
    and prints the name and the numbers within it.
    """
    # The predicted major product of the reaction is o-terphenyl.
    # Its systematic name is 1,2-diphenylbenzene.
    compound_name = "o-terphenyl"
    
    try:
        # Search for the compound by its common name in the PubChem database.
        compounds = pcp.get_compounds(compound_name, 'name')
        if compounds:
            # Get the first compound from the search results.
            product = compounds[0]
            
            # Retrieve the official IUPAC name.
            iupac_name = product.iupac_name
            
            print(f"The predicted product is known as '{compound_name}'.")
            print(f"Its official IUPAC name is: {iupac_name}")
            
            # Extract and print the numbers from the IUPAC name as requested.
            numbers = [char for char in iupac_name if char.isdigit()]
            print("The numbers in the IUPAC name are:")
            for num in numbers:
                print(num)
        else:
            print(f"Could not find the compound '{compound_name}' in the PubChem database.")
            print("The likely IUPAC name is 1,2-diphenylbenzene.")
            print("The numbers in the IUPAC name are:")
            print("1")
            print("2")

    except Exception as e:
        print(f"An error occurred. Please ensure 'pubchempy' is installed (`pip install pubchempy`).")
        # Provide a fallback answer if the library fails.
        print("The likely IUPAC name is 1,2-diphenylbenzene.")
        print("The numbers in the IUPAC name are:")
        print("1")
        print("2")

if __name__ == "__main__":
    find_iupac_name()