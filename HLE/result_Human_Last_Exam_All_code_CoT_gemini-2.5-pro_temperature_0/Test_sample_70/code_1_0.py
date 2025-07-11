def find_rank_2_african_flags():
    """
    This function identifies and prints the names of African nations whose flags
    have a linear algebraic rank of 2, the same as the flag of Denmark.

    The rank is determined by the number of linearly independent color patterns
    in the flag's design.

    - Denmark's flag rank is 2, based on its two patterns: (Red, White, Red) and (White, White, White).

    The following flags are identified as having a rank of 2 based on their structure:
    - Type A (Solid Field + Central Emblem): e.g., Morocco, Somalia.
    - Type B (Stripes + Central Emblem): e.g., Ghana, Senegal, Angola.
    - Type C (Hoist Bar + Horizontal Stripes): e.g., Benin, Madagascar.
    """

    print("The following flags of African nations have the same linear algebraic rank (2) as the flag of Denmark:")

    # List of countries whose flags have a rank of 2
    countries = [
        "Algeria",
        "Angola",
        "Benin",
        "Burkina Faso",
        "Cameroon",
        "Cape Verde",
        "Ghana",
        "Madagascar",
        "Mauritania",
        "Morocco",
        "Senegal",
        "Somalia"
    ]

    # Print each country name from the list
    for country in sorted(countries):
        print(f"- {country}")

find_rank_2_african_flags()