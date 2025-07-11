import collections

def solve_flag_rank_problem():
    """
    This script identifies African nations whose flags have the same
    maximal linear algebraic rank as the flag of Denmark.

    The maximal rank of a flag is equivalent to the number of distinct
    colors it contains. The flag of Denmark has 2 colors (red and white),
    so its rank is 2.

    The script finds all African flags that also have exactly 2 colors.
    """
    
    # A dictionary mapping African countries to a list of their flag colors.
    # Note: Variations in shades (e.g., light blue vs. blue) are treated as distinct.
    african_flags_colors = {
        "Algeria": ["Green", "White", "Red"],
        "Angola": ["Red", "Black", "Yellow"],
        "Benin": ["Green", "Yellow", "Red"],
        "Botswana": ["Light Blue", "Black", "White"],
        "Burkina Faso": ["Red", "Green", "Yellow"],
        "Burundi": ["Red", "Green", "White"],
        "Cabo Verde": ["Blue", "White", "Red", "Yellow"],
        "Cameroon": ["Green", "Red", "Yellow"],
        "Central African Republic": ["Blue", "White", "Green", "Yellow", "Red"],
        "Chad": ["Blue", "Yellow", "Red"],
        "Comoros": ["Yellow", "White", "Red", "Blue", "Green"],
        "Congo, Dem. Rep. of": ["Sky Blue", "Yellow", "Red"],
        "Congo, Rep. of the": ["Green", "Yellow", "Red"],
        "Cote d'Ivoire": ["Orange", "White", "Green"],
        "Djibouti": ["Light Blue", "Green", "White", "Red"],
        "Egypt": ["Red", "White", "Black", "Gold"],
        "Equatorial Guinea": ["Green", "White", "Red", "Blue"],
        "Eritrea": ["Green", "Red", "Blue", "Yellow"],
        "Eswatini": ["Blue", "Yellow", "Red", "Black", "White"],
        "Ethiopia": ["Green", "Yellow", "Red", "Blue"],
        "Gabon": ["Green", "Yellow", "Blue"],
        "Gambia, The": ["Red", "Blue", "Green", "White"],
        "Ghana": ["Red", "Yellow", "Green", "Black"],
        "Guinea": ["Red", "Yellow", "Green"],
        "Guinea-Bissau": ["Red", "Yellow", "Green", "Black"],
        "Kenya": ["Black", "Red", "Green", "White"],
        "Lesotho": ["Blue", "White", "Green", "Black"],
        "Liberia": ["Red", "White", "Blue"],
        "Libya": ["Red", "Black", "Green", "White"],
        "Madagascar": ["White", "Red", "Green"],
        "Malawi": ["Black", "Red", "Green"],
        "Mali": ["Green", "Yellow", "Red"],
        "Mauritania": ["Green", "Yellow", "Red"],
        "Mauritius": ["Red", "Blue", "Yellow", "Green"],
        "Morocco": ["Red", "Green"],
        "Mozambique": ["Green", "Black", "Yellow", "White", "Red"],
        "Namibia": ["Blue", "Red", "Green", "White", "Yellow"],
        "Niger": ["Orange", "White", "Green"],
        "Nigeria": ["Green", "White"],
        "Rwanda": ["Light Blue", "Yellow", "Green"],
        "Sao Tome and Principe": ["Green", "Yellow", "Red", "Black"],
        "Senegal": ["Green", "Yellow", "Red"],
        "Seychelles": ["Blue", "Yellow", "Red", "White", "Green"],
        "Sierra Leone": ["Green", "White", "Blue"],
        "Somalia": ["Light Blue", "White"],
        "South Africa": ["Red", "White", "Blue", "Green", "Black", "Yellow"],
        "South Sudan": ["Black", "Red", "Green", "White", "Blue", "Yellow"],
        "Sudan": ["Red", "White", "Black", "Green"],
        "Tanzania": ["Green", "Black", "Blue", "Yellow"],
        "Togo": ["Green", "Yellow", "Red", "White"],
        "Tunisia": ["Red", "White"],
        "Uganda": ["Black", "Yellow", "Red"],
        "Zambia": ["Green", "Red", "Black", "Orange"],
        "Zimbabwe": ["Green", "Yellow", "Red", "Black", "White"],
    }

    denmark_colors = ["Red", "White"]
    rank_denmark = len(set(denmark_colors))

    matching_countries = []
    for country, colors in african_flags_colors.items():
        # The rank is the number of unique colors
        rank_african_flag = len(set(colors))
        if rank_african_flag == rank_denmark:
            matching_countries.append(country)

    print(f"The maximal rank of the flag of Denmark is determined by its number of colors.")
    print(f"Rank(Denmark) = Number of colors in Danish flag = {rank_denmark}")
    print("\nFinding African flags with the same rank:")
    
    if not matching_countries:
        print("No African flags were found with the same rank as Denmark's flag.")
    else:
        for country in sorted(matching_countries):
            num_colors = len(set(african_flags_colors[country]))
            print(f"- {country}: Rank({country}) = Number of colors = {num_colors}")

    print("\nTherefore, the African nations with flags having the same rank as Denmark are:")
    print(", ".join(sorted(matching_countries)))

solve_flag_rank_problem()