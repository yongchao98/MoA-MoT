def find_matilda_female_colleagues():
    # Sarah's circle (including Sarah)
    sarah_males = 3
    sarah_females = 1  # This includes Matilda
    
    # John's circle (including John)
    john_males = 3
    john_females = 2  # This includes Matilda
    
    # We know Matilda is in both circles
    # Let's count her female colleagues:
    
    # In Sarah's circle:
    # - Sarah herself is a female colleague
    # - No other females (as the 1 female mentioned is Matilda)
    females_from_sarah = 1  # Just Sarah
    
    # In John's circle:
    # - The other female besides Matilda
    females_from_john = 1  # The other female in John's circle
    
    total_female_colleagues = females_from_sarah + females_from_john
    
    print(f"Matilda has {total_female_colleagues} female colleagues")
    return total_female_colleagues

find_matilda_female_colleagues()