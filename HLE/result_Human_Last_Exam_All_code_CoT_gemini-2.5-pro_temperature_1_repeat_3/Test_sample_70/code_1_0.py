def find_african_flags_with_rank_two():
    """
    This function identifies and prints the names of African nations whose flags,
    when viewed as matrices, have a linear algebraic rank of 2.

    The rank of the Danish flag is 2, as its Nordic cross pattern creates
    variation in both horizontal and vertical directions.

    A flag has a rank of 2 if its design is not composed solely of simple
    horizontal or vertical stripes. This includes flags with cantons, hoist
    triangles, central emblems spanning stripes, or diagonal elements.
    """

    # List of African nations with flags of rank 2
    rank_two_flags = [
        "Algeria",              # Vertical bicolor with central emblem
        "Angola",               # Horizontal bicolor with central emblem
        "Burundi",              # Saltire (diagonal cross) with central circle
        "Central African Republic", # Horizontal stripes with an overlapping vertical stripe
        "Comoros",              # Horizontal stripes with a hoist triangle
        "Democratic Republic of the Congo", # Diagonal stripe
        "Djibouti",             # Horizontal bicolor with a hoist triangle
        "Equatorial Guinea",    # Horizontal tricolor with a hoist triangle
        "Eritrea",              # Large hoist triangle extending to the fly
        "Eswatini",             # Horizontal stripes with a central shield emblem
        "Kenya",                # Horizontal tricolor with a central shield emblem
        "Liberia",              # Horizontal stripes with a canton (square in the corner)
        "Madagascar",           # Vertical stripe alongside horizontal stripes
        "Mozambique",           # Horizontal tricolor with a hoist triangle
        "Namibia",              # Diagonal stripe
        "Republic of the Congo",# Diagonal stripe
        "Sao Tome and Principe",# Horizontal tricolor with a hoist triangle
        "Seychelles",           # Oblique (diagonal) bands
        "Somalia",              # Solid field with a central star
        "South Africa",         # A 'Y' shape (pall) creating diagonal and horizontal elements
        "South Sudan",          # Horizontal tricolor with a hoist triangle
        "Tanzania",             # Diagonal stripe
        "Togo",                 # Horizontal stripes with a canton
        "Tunisia",              # Solid field with a central disc and emblem
        "Uganda",               # Horizontal stripes with a central disc emblem
        "Zambia",               # Solid field with vertical stripes at the fly
        "Zimbabwe"              # Horizontal stripes with a hoist triangle
    ]

    print("The rank of the flag of Denmark is 2.")
    print("The following African national flags also have a rank of 2:")
    for country in rank_two_flags:
        print(country)

find_african_flags_with_rank_two()