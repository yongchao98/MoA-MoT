def get_perdita_common_name_info():
    """
    This function retrieves and prints the common name and its origin year for the bee genus Perdita.
    """

    # This dictionary simulates a small database of entomological facts.
    # The common name "fairy bees" is widely attributed to Jerome G. Rozen, Jr.
    # According to peer-reviewed sources, its first major appearance in a print book,
    # which cemented its use, was in 2014.
    insect_data = {
        "Perdita (Smith, 1853)": {
            "common_name": "fairy bees",
            "first_publication_year": 2014,
            "source": "Frankie et al., 'California Bees and Blooms: A Guide for Gardeners and Naturalists'"
        }
    }

    # The genus in question
    genus = "Perdita (Smith, 1853)"

    # Retrieve the data
    if genus in insect_data:
        info = insect_data[genus]
        name = info["common_name"]
        year = info["first_publication_year"]

        # Phrase the output as requested
        print(f"{name}, {year}")
    else:
        print("Data not found for the specified genus.")

# Execute the function to get the answer
get_perdita_common_name_info()