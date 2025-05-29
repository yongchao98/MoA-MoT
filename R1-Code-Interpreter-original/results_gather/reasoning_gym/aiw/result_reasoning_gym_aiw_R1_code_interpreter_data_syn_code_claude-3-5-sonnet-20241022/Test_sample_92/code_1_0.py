def count_matilda_female_colleagues():
    # From Sarah's circle
    sarah_circle_females = 1  # Only Sarah, as Matilda is the other female
    
    # From James' circle
    james_circle_females = 3  # Total females in James' circle
    james_female_colleagues = james_circle_females - 1  # Subtract Matilda herself
    
    # Total female colleagues of Matilda
    total = sarah_circle_females + james_female_colleagues
    
    print(total)

count_matilda_female_colleagues()