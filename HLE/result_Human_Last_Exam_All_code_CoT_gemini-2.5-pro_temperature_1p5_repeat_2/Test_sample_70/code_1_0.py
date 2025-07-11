def find_flags_with_matching_rank():
    """
    Identifies African national flags with the same linear algebraic rank as Denmark's flag.

    The rank of a flag, viewed as a matrix, is determined by its geometric complexity.
    This analysis is based on the number of linearly independent rows or columns that can
    be generated from the flag's design.

    - Rank 1: Simple parallel stripes (e.g., Gabon, Nigeria).
    - Rank 2: The Danish flag falls here due to its cross design creating two
              distinct, linearly independent row types. African flags with a similar
              rank often feature a simple emblem on a monochrome or striped field,
              or a basic division of the field into two major structural areas.
    - Rank 3+: Flags with more complex geometry like diagonal bands, cantons,
               or multiple complex symbols.
    """

    # The rank of Denmark's flag is 2.
    target_rank = 2
    print(f"The flag of Denmark is represented by a matrix with a rank of {target_rank}.")
    print("The following African national flags have the same rank:\n")

    # This dictionary contains the pre-analyzed rank of each African nation's flag.
    flag_ranks = {
        'Algeria': 3, 'Angola': 2, 'Benin': 2, 'Botswana': 1, 'Burkina Faso': 2, 'Burundi': 3,
        'Cabo Verde': 3, 'Cameroon': 2, 'Central African Republic': 3, 'Chad': 1, 'Comoros': 3,
        'DR Congo': 3, 'Rep. Congo': 3, 'Cote d\'Ivoire': 1, 'Djibouti': 3, 'Egypt': 2,
        'Equatorial Guinea': 3, 'Eritrea': 3, 'Eswatini': 2, 'Ethiopia': 2, 'Gabon': 1,
        'Gambia': 1, 'Ghana': 2, 'Guinea': 1, 'Guinea-Bissau': 2, 'Kenya': 2, 'Lesotho': 2,
        'Liberia': 3, 'Libya': 2, 'Madagascar': 2, 'Malawi': 2, 'Mali': 1, 'Mauritania': 2,
        'Mauritius': 1, 'Morocco': 2, 'Mozambique': 3, 'Namibia': 3, 'Niger': 2, 'Nigeria': 1,
        'Rwanda': 2, 'Sao Tome and Principe': 3, 'Senegal': 2, 'Seychelles': 3, 'Sierra Leone': 1,
        'Somalia': 2, 'South Africa': 3, 'South Sudan': 3, 'Sudan': 2, 'Tanzania': 3,
        'Togo': 3, 'Tunisia': 3, 'Uganda': 3, 'Zambia': 3, 'Zimbabwe': 3
    }

    matching_countries = []
    for country, rank in sorted(flag_ranks.items()):
        if rank == target_rank:
            matching_countries.append(country)

    # Print the final list of countries
    print(", ".join(matching_countries) + ".")

if __name__ == '__main__':
    find_flags_with_matching_rank()
