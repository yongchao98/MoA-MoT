def analyze_colleagues():
    # Linda's circle (including Linda)
    linda_males = 4
    linda_females = 2  # excluding Linda
    linda_total = linda_males + linda_females + 1  # +1 for Linda

    # David's circle (including David)
    david_males = 4
    david_females = 4  # excluding David
    david_total = david_males + david_females + 1  # +1 for David

    # We know Matilda is in both circles
    # Let's find Matilda's female colleagues:
    
    # In Linda's circle:
    # - Linda herself
    # - The other female colleague besides Matilda
    matilda_females_from_linda = 2  # Linda + 1 other female
    
    # In David's circle:
    # - The 4 females (including Matilda)
    # So there are 3 other females in David's circle
    matilda_females_from_david = 3  # excluding Matilda
    
    # Total unique female colleagues of Matilda
    total_matilda_females = matilda_females_from_linda + matilda_females_from_david
    
    print(f"Matilda has {total_matilda_females} female colleagues")

analyze_colleagues()