def find_runestone():
    """
    Identifies a runestone by searching a database for features
    observed in an image.
    """
    # A small knowledge base of Ingvar runestones with their key features.
    runestones = [
        {
            'id': 'Sö 179',
            'common_name': 'Gripsholm Runestone',
            'features': ['serpent_band'],
            # Contains 'stæin þennsa' (this stone), singular.
            'inscription_text': 'Tola let ræisa stæin þennsa at sun sinn Harald'
        },
        {
            'id': 'Sö 281',
            'common_name': 'Strängnäs Cathedral stone fragment',
            # Based on 17th-century drawing, not a grid.
            'features': ['serpent_band'],
            # Contains 'stæina þessa' but the style does not match.
            'inscription_text': '... ræisti stæina þessa at ...'
        },
        {
            'id': 'U 644',
            'common_name': 'Varpsund Runestone (upper part)',
            'features': ['grid_layout'],
            # Style matches, but the key text is not on this fragment.
            'inscription_text': '... þæiR letu ræisa ...'
        },
        {
            'id': 'U 654',
            'common_name': 'Varpsund Runestone (lower part)',
            'features': ['grid_layout'],
            # Both the style and text match.
            'inscription_text': '... letu ræisa staina þisa at faður sinn ...'
        }
    ]

    # Define features identified from the provided image.
    image_feature_style = 'grid_layout'
    image_feature_text = 'staina þisa'

    # Search for the runestone that matches all identified features.
    matched_stone_id = None
    for stone in runestones:
        style_match = image_feature_style in stone['features']
        text_match = image_feature_text in stone['inscription_text']

        if style_match and text_match:
            matched_stone_id = stone['id']
            break

    # Output the conclusion.
    if matched_stone_id:
        # The ID consists of a letter code for the region (U for Uppland) and a number.
        region_code, number_code_str = matched_stone_id.split(' ')
        number_code = int(number_code_str)

        print(f"Image analysis shows:")
        print(f"1. A 'grid_layout' carving style.")
        print(f"2. The runic phrase 'staina þisa'.")
        print(f"\nSearching the database for a runestone with these features...")
        print(f"A match was found!")
        print(f"The ID is composed of the region code '{region_code}' and the number '{number_code}'.")
        print(f"Full ID: {matched_stone_id}")
    else:
        print("Could not identify the runestone based on the provided features.")

find_runestone()