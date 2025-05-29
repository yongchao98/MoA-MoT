def analyze_colleagues():
    # Jennifer's circle (including Jennifer)
    jennifer_males = 2
    jennifer_females = 4  # excluding Jennifer
    
    # John's circle (including John)
    john_males = 3
    john_females = 6  # excluding John
    
    # Since Matilda is in both circles and is female
    # Jennifer's actual female colleagues = total - 1 (Matilda)
    jennifer_unique_females = jennifer_females - 1
    
    # For Matilda:
    # From Jennifer's circle: all females (except herself) + all males + Jennifer
    matilda_from_jennifer = jennifer_unique_females + jennifer_males + 1
    
    # From John's circle: all females (except herself) + all males + John
    matilda_from_john = (john_females - 1) + john_males + 1
    
    # Total unique colleagues for Matilda
    total_matilda_colleagues = matilda_from_jennifer + matilda_from_john
    
    # Calculate female colleagues
    # From Jennifer's circle: jennifer_unique_females + Jennifer = jennifer_females
    # From John's circle: john_females - 1 (excluding Matilda)
    matilda_female_colleagues = jennifer_unique_females + 1 + (john_females - 1)
    
    print(f"Matilda has {matilda_female_colleagues} female colleagues")

analyze_colleagues()