def analyze_colleagues():
    # Sarah's circle (including Sarah)
    sarah_males = 5
    sarah_females = 5  # excluding Sarah herself
    
    # John's circle (including John)
    john_males = 2     # excluding John himself
    john_females = 4
    
    # We know Matilda is in both circles
    # Let's count Matilda's female colleagues:
    
    # From Sarah's circle:
    # All females in Sarah's circle (excluding Matilda herself)
    females_from_sarah = sarah_females - 1  # -1 for Matilda
    
    # From John's circle:
    # All females in John's circle (excluding Matilda herself)
    females_from_john = john_females - 1    # -1 for Matilda
    
    # Since people in John's circle don't have other colleagues
    # (except Matilda), the females in John's circle are a subset
    # of the females in Sarah's circle
    
    # Total female colleagues of Matilda:
    total_females = females_from_sarah
    
    print(f"Matilda has {total_females} female colleagues")

analyze_colleagues()