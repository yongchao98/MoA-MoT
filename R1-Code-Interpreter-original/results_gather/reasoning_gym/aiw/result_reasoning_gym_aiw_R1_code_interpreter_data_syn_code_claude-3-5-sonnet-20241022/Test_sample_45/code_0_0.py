def find_matilda_female_colleagues():
    # Barbara's circle (excluding Barbara)
    barbara_males = 2
    barbara_females = 2  # including Matilda
    
    # Michael's circle (excluding Michael)
    michael_males = 1
    michael_females = 3  # including Matilda
    
    # Since Matilda is in both circles, let's count her unique female colleagues
    
    # In Barbara's circle (excluding Matilda herself)
    barbara_circle_females = barbara_females - 1  # one is Matilda
    
    # In Michael's circle (excluding Matilda herself)
    michael_circle_females = michael_females - 1  # one is Matilda
    
    # Total unique female colleagues of Matilda
    total_female_colleagues = barbara_circle_females + michael_circle_females
    
    print(f"Matilda's female colleagues: {total_female_colleagues}")

find_matilda_female_colleagues()