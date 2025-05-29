def analyze_colleagues():
    # Patricia's circle (excluding Patricia)
    patricia_males = 6
    patricia_females = 6
    
    # Charles's circle (excluding Charles)
    charles_males = 6
    charles_females = 5
    
    # We know Matilda is in both circles and is female
    # So she's counted in both patricia_females and charles_females
    
    # In Patricia's circle, Matilda has:
    matilda_colleagues_from_patricia = patricia_males + (patricia_females - 1)  # -1 because we don't count Matilda herself
    
    # In Charles's circle, Matilda has:
    matilda_colleagues_from_charles = charles_males + (charles_females - 1)  # -1 because we don't count Matilda herself
    
    # Since Matilda is in both circles, we need to count unique female colleagues
    # From Patricia's circle: patricia_females - 1 (excluding herself)
    # From Charles's circle: charles_females - 1 (excluding herself)
    total_female_colleagues = patricia_females - 1  # Start with Patricia's circle
    
    # Add any unique females from Charles's circle
    # We know there's overlap of exactly 1 female (Matilda herself)
    # So we add remaining females from Charles's circle
    total_female_colleagues += (charles_females - 1)
    
    print(total_female_colleagues)

analyze_colleagues()