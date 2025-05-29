def verify_matilda_colleagues():
    # Known facts
    charles_female_colleagues = 5  # excluding Charles
    barbara_female_colleagues = 1  # excluding Barbara
    
    # Calculate Matilda's female colleagues
    matilda_females = 0
    
    # From Charles' circle
    matilda_females += (charles_female_colleagues - 1)  # -1 because Matilda is counted in this number
    
    # From Barbara's circle
    matilda_females += 1  # Barbara herself
    
    # No need to add barbara_female_colleagues because Matilda is that one female colleague
    
    print(f"Detailed count:")
    print(f"Females from Charles' circle (excluding Matilda): {charles_female_colleagues - 1}")
    print(f"Barbara herself: 1")
    print(f"Total female colleagues of Matilda: {matilda_females}")

verify_matilda_colleagues()