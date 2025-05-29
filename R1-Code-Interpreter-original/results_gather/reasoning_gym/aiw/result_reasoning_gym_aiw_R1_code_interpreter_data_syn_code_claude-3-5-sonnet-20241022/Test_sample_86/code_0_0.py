def analyze_colleagues():
    # Jennifer's circle (including Jennifer)
    jennifer_males = 1
    jennifer_females = 6  # other females excluding Jennifer
    jennifer_total = jennifer_males + jennifer_females + 1  # +1 for Jennifer herself
    
    # William's circle (including William)
    william_males = 3
    william_females = 3  # other females excluding William
    william_total = william_males + william_females + 1  # +1 for William himself
    
    # We know Matilda is in both circles and is female
    # So she's counted in both jennifer_females and william_females
    
    # Matilda's female colleagues would be:
    # From Jennifer's circle: all females except herself = jennifer_females
    # From William's circle: all females except herself = william_females
    # Need to avoid double counting
    
    matilda_female_colleagues = jennifer_females + william_females - 1  # -1 because Matilda is counted twice
    
    print(f"Matilda's female colleagues: {matilda_female_colleagues}")

analyze_colleagues()