def check_paintings():
    """
    This function analyzes a predefined list of paintings to determine
    which ones were created by a French painter before the year 1900.
    The data is based on art historical research for the provided images.
    """
    paintings_data = [
        {'label': 'A', 'nationality': 'French', 'year': 1888},
        {'label': 'B', 'nationality': 'French', 'year': 1898},
        {'label': 'C', 'nationality': 'French', 'year': 1892},
        {'label': 'D', 'nationality': 'French', 'year': 1890},
        {'label': 'E', 'nationality': 'French', 'year': 1903},
        {'label': 'F', 'nationality': 'N/A', 'year': 2022}
    ]

    result = []
    for painting in paintings_data:
        is_french_before_1900 = painting['nationality'] == 'French' and painting['year'] < 1900
        if is_french_before_1900:
            result.append('Y')
        else:
            result.append('N')

    print("".join(result))

check_paintings()