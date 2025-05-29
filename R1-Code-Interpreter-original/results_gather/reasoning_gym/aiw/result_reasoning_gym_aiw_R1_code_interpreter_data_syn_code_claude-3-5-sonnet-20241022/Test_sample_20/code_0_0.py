def analyze_colleagues():
    # Mary's circle (including Mary)
    mary_males = 4
    mary_females = 3  # excluding Mary herself
    
    # Thomas's circle (including Thomas)
    thomas_males = 4
    thomas_females = 1  # excluding Thomas himself
    
    # We know Matilda is in both circles and is female
    # Let's count Matilda's female colleagues:
    
    # In Mary's circle:
    # - All other females in Mary's circle (excluding Matilda)
    females_from_mary = mary_females - 1  # -1 because Matilda is counted in mary_females
    
    # In Thomas's circle:
    # - No additional females because Thomas has only 1 female colleague (Matilda)
    females_from_thomas = 0  # All females in Thomas's circle is just Matilda
    
    total_female_colleagues = females_from_mary + females_from_thomas
    
    print(f"Matilda has {total_female_colleagues} female colleagues")

analyze_colleagues()