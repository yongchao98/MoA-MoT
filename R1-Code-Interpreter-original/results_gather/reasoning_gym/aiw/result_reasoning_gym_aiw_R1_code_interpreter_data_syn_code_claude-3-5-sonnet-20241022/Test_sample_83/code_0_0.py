def analyze_colleagues():
    # Linda's circle (including Linda)
    linda_males = 5
    linda_females = 1  # This includes Matilda
    linda_total = linda_males + linda_females + 1  # +1 for Linda herself
    
    # John's circle (including John)
    john_males = 3
    john_females = 5  # This includes Matilda
    john_total = john_males + john_females + 1  # +1 for John himself
    
    # Since Matilda is in both circles, let's count her female colleagues
    
    # In Linda's circle:
    # - Linda herself is a female colleague
    matilda_females_from_linda = 1  # Just Linda
    
    # In John's circle:
    # - All females except Matilda herself
    matilda_females_from_john = john_females - 1
    
    # Total female colleagues of Matilda
    total_matilda_females = matilda_females_from_linda + matilda_females_from_john
    
    print(f"Matilda has {total_matilda_females} female colleagues")
    return total_matilda_females

analyze_colleagues()