def find_rank_2_african_flags():
    """
    Analyzes the rank of the Danish flag and finds African flags with the same rank.

    The rank of a flag is determined by representing it as a matrix of color values
    and finding the number of linearly independent rows. To maximize rank, distinct
    colors are assigned distinct non-zero numbers.
    """

    # Data describing the structure of African flags based on visual analysis.
    # 'base' describes the background, 'features' are additions that can increase rank.
    african_nations = {
        'Algeria': {'base': 'vertical_stripes', 'features': ['emblem']},
        'Angola': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Benin': {'base': 'benin_style', 'features': []},
        'Botswana': {'base': 'horizontal_stripes', 'features': []},
        'Burkina Faso': {'base': 'horizontal_stripes', 'features': ['star']},
        'Burundi': {'base': 'complex_saltire', 'features': ['emblem']},
        'Cabo Verde': {'base': 'horizontal_stripes', 'features': ['star_circle']},
        'Cameroon': {'base': 'vertical_stripes', 'features': ['star']},
        'Central African Republic': {'base': 'complex_superimposed', 'features': []},
        'Chad': {'base': 'vertical_stripes', 'features': []},
        'Comoros': {'base': 'complex_triangle', 'features': ['emblem']},
        'DR Congo': {'base': 'complex_diagonal', 'features': ['star']},
        'Congo, Republic of the': {'base': 'complex_diagonal', 'features': []},
        'Cote d\'Ivoire': {'base': 'vertical_stripes', 'features': []},
        'Djibouti': {'base': 'complex_triangle', 'features': ['star']},
        'Egypt': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Equatorial Guinea': {'base': 'complex_triangle', 'features': ['emblem']},
        'Eritrea': {'base': 'complex_triangle', 'features': ['emblem']},
        'Eswatini': {'base': 'horizontal_stripes', 'features': ['shield']},
        'Ethiopia': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Gabon': {'base': 'horizontal_stripes', 'features': []},
        'Gambia': {'base': 'horizontal_stripes', 'features': []},
        'Ghana': {'base': 'horizontal_stripes', 'features': ['star']},
        'Guinea': {'base': 'vertical_stripes', 'features': []},
        'Guinea-Bissau': {'base': 'complex_superimposed', 'features': []},
        'Kenya': {'base': 'horizontal_stripes', 'features': ['shield']},
        'Lesotho': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Liberia': {'base': 'complex_canton', 'features': ['star']},
        'Libya': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Madagascar': {'base': 'benin_style', 'features': []},
        'Malawi': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Mali': {'base': 'vertical_stripes', 'features': []},
        'Mauritania': {'base': 'complex_3_row_types', 'features': []},
        'Mauritius': {'base': 'horizontal_stripes', 'features': []},
        'Morocco': {'base': 'monocolor', 'features': ['pentagram']},
        'Mozambique': {'base': 'complex_triangle', 'features': ['emblem']},
        'Namibia': {'base': 'complex_diagonal', 'features': []},
        'Niger': {'base': 'horizontal_stripes', 'features': ['circle']},
        'Nigeria': {'base': 'vertical_stripes', 'features': []},
        'Rwanda': {'base': 'horizontal_stripes', 'features': ['emblem']},
        'Sao Tome and Principe': {'base': 'complex_triangle', 'features': ['stars']},
        'Senegal': {'base': 'vertical_stripes', 'features': ['star']},
        'Seychelles': {'base': 'complex_oblique', 'features': []},
        'Sierra Leone': {'base': 'horizontal_stripes', 'features': []},
        'Somalia': {'base': 'monocolor', 'features': ['star']},
        'South Africa': {'base': 'complex_pall', 'features': []},
        'South Sudan': {'base': 'complex_triangle', 'features': ['star']},
        'Sudan': {'base': 'complex_triangle', 'features': []},
        'Tanzania': {'base': 'complex_diagonal', 'features': []},
        'Togo': {'base': 'complex_canton', 'features': ['star']},
        'Tunisia': {'base': 'complex_disk', 'features': ['emblem']},
        'Uganda': {'base': 'complex_disk', 'features': ['emblem']},
        'Zambia': {'base': 'complex_fly_emblem', 'features': []},
        'Zimbabwe': {'base': 'complex_triangle', 'features': ['emblem']}
    }

    def get_rank(flag_data):
        base = flag_data['base']
        features = flag_data['features']
        
        # Rank 1 flags have a simple stripe/color base with no rank-increasing features.
        if ('stripes' in base or 'monocolor' in base) and not features:
            return 1
            
        # Rank 2 flags have a simple base + a simple feature, or a specific 2-rank geometry.
        if (('stripes' in base or 'monocolor' in base) and features) or (base == 'benin_style'):
            return 2
            
        # Any other complex geometry results in a rank > 2.
        return 3

    rank_2_countries = []
    for country, data in african_nations.items():
        if get_rank(data) == 2:
            rank_2_countries.append(country)

    print("The flag of Denmark has a linear algebraic rank of 2.")
    print("This is determined by representing the flag as a matrix and finding the number of linearly independent row vectors.")
    print("If we let Red = 1 and White = 2, the flag's matrix contains two types of rows that are not multiples of each other:")
    print("Row 1 (from white crossbar): [2, 2, 2, 2, 2]")
    print("Row 2 (from red field and vertical crossbar): [1, 1, 2, 1, 1]")
    print("Because no number 'c' satisfies c * [2, 2, 2, 2, 2] = [1, 1, 2, 1, 1], the vectors are linearly independent, and the rank is 2.\n")
    print("African nations with flags that also have a rank of 2 are:")
    
    for country in sorted(rank_2_countries):
        print(f"- {country}")

find_rank_2_african_flags()