def find_rank_two_flags():
    """
    Prints the list of African nations whose flags have the same linear algebraic rank (2) as the flag of Denmark.

    The determination of rank is based on the flag's geometric structure:
    - Rank 1: Simple parallel stripes (all horizontal or all vertical).
    - Rank 2: A cross, stripes with an emblem/canton, a field with an emblem, a hoist-side element, etc.
    - Rank 3+: More complex designs like Y-shapes, saltires, or multiple diagonal elements.
    """

    # List of African nations with flags of rank 2
    rank_two_nations = [
        "Algeria",
        "Angola",
        "Burkina Faso",
        "Cameroon",
        "Cape Verde",
        "Central African Republic",
        "Comoros",
        "Djibouti",
        "Egypt",
        "Equatorial Guinea",
        "Eswatini",
        "Ethiopia",
        "Ghana",
        "Kenya",
        "Lesotho",
        "Liberia",
        "Libya",
        "Madagascar",
        "Malawi",
        "Mauritania",
        "Morocco",
        "Mozambique",
        "Niger",
        "Rwanda",
        "Sao Tome and Principe",
        "Senegal",
        "Somalia",
        "South Sudan",
        "Sudan",
        "Togo",
        "Tunisia",
        "Uganda",
        "Zambia",
        "Zimbabwe"
    ]

    print("The flags of the following African nations have the same linear algebraic rank as the flag of Denmark (Rank 2):")
    for nation in sorted(rank_two_nations):
        print(f"- {nation}")

if __name__ == "__main__":
    find_rank_two_flags()
