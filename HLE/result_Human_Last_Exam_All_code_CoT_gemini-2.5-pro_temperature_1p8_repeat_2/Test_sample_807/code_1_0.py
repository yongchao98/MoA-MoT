def find_crayon_color():
    """
    Finds the Crayola crayon color that contains a specific chemical.
    """
    # The chemical in question is commonly known as Rhodamine B.
    chemical_data = {
        "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride": {
            "common_name": "Rhodamine B",
            "crayon_color": "Screamin' Green"
        }
    }

    query_chemical = "9-(2-Carboxyphenyl)-6-(diethylamino)-N,N-diethyl-3H-xanthen-3-iminium chloride"

    if query_chemical in chemical_data:
        info = chemical_data[query_chemical]
        color = info["crayon_color"]
        common_name = info["common_name"]

        print(f"The chemical provided is: {query_chemical}")
        print(f"Its common name is: {common_name}")
        print(f"\nThis pigment is used to create the Crayola crayon color: {color}")
        print("\nInterestingly, this pink fluorescent dye is mixed with a yellow pigment to produce a color that appears green in normal light but fluoresces orange-red under UV light.")
    else:
        print("Crayon color not found for the specified chemical.")

find_crayon_color()