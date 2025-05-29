def analyze_colleagues():
    # Linda's circle
    linda_males = 5
    linda_females = 3
    
    # Richard's circle
    richard_males = 2
    richard_females = 4
    
    # Since Matilda is in both circles:
    # She is one of Linda's female colleagues
    # She is one of Richard's female colleagues
    
    # Matilda's female colleagues would be:
    # All females in Linda's circle (except herself) +
    # All females in Richard's circle (except herself)
    # But we need to avoid double counting
    
    # In Linda's circle:
    matilda_females_from_linda = linda_females - 1  # excluding Matilda herself
    
    # In Richard's circle:
    matilda_females_from_richard = richard_females - 1  # excluding Matilda herself
    
    # Since Richard's circle people don't have other colleagues (except Matilda),
    # these must be different groups of people
    total_female_colleagues = matilda_females_from_linda + matilda_females_from_richard
    
    print(total_female_colleagues)

analyze_colleagues()