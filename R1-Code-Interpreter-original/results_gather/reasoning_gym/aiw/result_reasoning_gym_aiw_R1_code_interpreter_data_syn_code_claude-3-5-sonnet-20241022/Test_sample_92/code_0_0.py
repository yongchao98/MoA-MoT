def analyze_colleagues():
    # Sarah's circle (including Sarah)
    sarah_males = 2
    sarah_females = 1  # excluding Sarah herself
    sarah_total = sarah_males + sarah_females + 1  # +1 for Sarah

    # James' circle (including James)
    james_males = 4
    james_females = 3  # excluding James himself
    james_total = james_males + james_females + 1  # +1 for James

    # We know Matilda is in both circles
    # Let's find Matilda's female colleagues:
    
    # In Sarah's circle:
    # - Sarah herself (1 female)
    # - The other female colleague mentioned
    sarah_female_colleagues = 1 + 1  # Sarah + 1 other female
    
    # In James' circle:
    # - The 3 females mentioned (including Matilda herself)
    james_female_colleagues = 3  # Total females in James' circle
    
    # Since Matilda is counted in both circles
    # And she doesn't count herself as a colleague
    total_female_colleagues = (sarah_female_colleagues - 1) + (james_female_colleagues - 1)
    
    print(total_female_colleagues)

analyze_colleagues()