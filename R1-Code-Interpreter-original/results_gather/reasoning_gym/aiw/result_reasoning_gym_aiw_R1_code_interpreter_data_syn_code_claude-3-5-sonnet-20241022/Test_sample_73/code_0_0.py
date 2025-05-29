def analyze_colleagues():
    # Barbara's circle (including Barbara)
    barbara_males = 1
    barbara_females = 1  # excluding Barbara
    barbara_total = barbara_males + barbara_females + 1  # +1 for Barbara herself
    
    # Charles's circle (including Charles)
    charles_males = 3
    charles_females = 5  # excluding Charles
    charles_total = charles_males + charles_females + 1  # +1 for Charles himself
    
    # We know Matilda is in both circles
    # She is one of the females in both counts
    
    # Matilda's female colleagues would be:
    # - All females in Barbara's circle (except herself) = barbara_females
    # - All females in Charles's circle (except herself) = charles_females
    # - Plus Barbara if she's female (which she is)
    # But we need to avoid double counting
    
    matilda_female_colleagues = charles_females + barbara_females
    
    # If Matilda is counted in both circles, we need to subtract 1
    matilda_female_colleagues -= 1  # remove double counting of Matilda
    
    # Add Barbara
    matilda_female_colleagues += 1
    
    print(f"Matilda has {matilda_female_colleagues} female colleagues")

analyze_colleagues()