import re

def find_crayon_color():
    """
    Identifies the Crayola crayon color containing a specific chemical.
    The information is based on known pigment compositions.
    """
    # The chemical name provided by the user
    chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # This chemical is commonly known as Rhodamine B, a dye used for a specific fluorescent color.
    # We can store this known information in a dictionary to simulate a database lookup.
    pigment_database = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Red Violet"
    }

    # Retrieve the crayon color from our 'database'
    crayon_color = pigment_database.get(chemical_name, "Unknown")

    # Extract all numbers from the chemical name string as requested
    numbers = re.findall(r'\d+', chemical_name)

    # Print the results
    print(f"The chemical is: {chemical_name}")
    print(f"The numbers in the chemical name are: {', '.join(numbers)}")
    print("-" * 20)
    print(f"The Crayola crayon color containing this chemical is: {crayon_color}")

# Execute the function to find and print the answer
find_crayon_color()