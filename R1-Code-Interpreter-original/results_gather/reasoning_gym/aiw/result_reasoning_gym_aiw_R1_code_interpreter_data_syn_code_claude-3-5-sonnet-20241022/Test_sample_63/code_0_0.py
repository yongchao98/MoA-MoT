def analyze_colleagues():
    # Sarah's circle (including Sarah)
    sarah_males = 2
    sarah_females = 4  # excluding Sarah herself
    
    # Thomas's circle (including Thomas)
    thomas_males = 1
    thomas_females = 2  # excluding Thomas himself
    
    # We know Matilda is in both circles and is female
    # So she's one of Sarah's 4 female colleagues
    # And one of Thomas's 2 female colleagues
    
    # Matilda's female colleagues would be:
    # - All females in Sarah's circle (except herself) = sarah_females - 1
    # - All females in Thomas's circle (except herself) = thomas_females - 1
    # But we need to avoid double counting
    
    matilda_female_colleagues = (sarah_females - 1) + 0  # We don't add Thomas's females because they don't have outside connections except Matilda
    
    print(f"Matilda has {matilda_female_colleagues} female colleagues")

analyze_colleagues()