def analyze_colleagues():
    # Known facts
    susan_males = 2
    susan_females = 4
    charles_males = 4
    charles_females = 6
    
    # Since Matilda is female and in both circles
    # In Susan's circle:
    # - Matilda is one of the 4 females
    # - The other 3 females and 2 males are unique to Susan's circle
    
    # In Charles's circle:
    # - Matilda is one of the 6 females
    # - The other 5 females and 4 males are unique to Charles's circle
    
    # Matilda's colleagues:
    # From Susan's circle: 2 males + 3 females (excluding herself)
    # From Charles's circle: 4 males + 5 females (excluding herself)
    # Total unique females = 3 (from Susan) + 5 (from Charles) = 8
    
    matilda_female_colleagues = 8
    print(f"Matilda has {matilda_female_colleagues} female colleagues")
    
analyze_colleagues()