def calculate_matilda_female_colleagues():
    # Susan's circle (excluding Matilda)
    susan_males = 5
    susan_females = 4  # 5-1 because Matilda is one of them
    
    # James's circle (excluding Matilda)
    james_males = 2
    james_females = 5  # 6-1 because Matilda is one of them
    
    # Total female colleagues of Matilda:
    # - All females in Susan's circle (except herself)
    # - All females in James's circle (except herself)
    # - But we need to avoid double counting if there are common people
    # Since Matilda is the only connection between circles,
    # there are no other common people
    
    total_female_colleagues = susan_females + james_females
    
    print(f"Matilda's female colleagues: {total_female_colleagues}")
    return total_female_colleagues

calculate_matilda_female_colleagues()