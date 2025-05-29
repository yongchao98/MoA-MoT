def analyze_colleagues():
    # Patricia's circle (including Patricia)
    patricia_males = 1
    patricia_females = 5  # excluding Patricia
    patricia_total = patricia_males + patricia_females + 1  # +1 for Patricia

    # Robert's circle (including Robert)
    robert_males = 3  # excluding Robert
    robert_females = 3
    robert_total = robert_males + robert_females + 1  # +1 for Robert

    # We know Matilda is in both circles and she's female
    # So she's counted in both patricia_females and robert_females

    # Matilda's female colleagues would be:
    # From Patricia's circle: all females except herself = patricia_females - 1
    # Plus Patricia herself = +1
    matilda_females_from_patricia = patricia_females

    # From Robert's circle: all females except herself = robert_females - 1
    # Plus any additional females not counted yet
    matilda_females_from_robert = robert_females - 1  # -1 because Matilda is counted in robert_females

    # Total female colleagues of Matilda (removing duplicates)
    # We need to combine both circles but be careful not to double count
    total_matilda_females = matilda_females_from_patricia  # All females from Patricia's circle
    
    print(f"Matilda's female colleagues: {total_matilda_females}")

analyze_colleagues()