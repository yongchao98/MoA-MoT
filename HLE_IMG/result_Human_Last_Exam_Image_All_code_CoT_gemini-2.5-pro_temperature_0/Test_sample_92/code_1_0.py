def solve_art_puzzle():
    """
    Determines which paintings were created by a French painter before 1900.
    The logic is based on art historical facts for each image.
    """
    
    # Data for each painting: (ID, Artist Nationality, Year of Creation)
    paintings = {
        'A': ('Swiss', 1892),   # Félix Vallotton, "Summer Evening, Bathers"
        'B': ('French', 1879),  # Edgar Degas, "At the Louvre: The Etruscan Sarcophagus"
        'C': ('French', 1891),  # Henri-Edmond Cross, "The Olive Trees"
        'D': ('French', 1890),  # Odilon Redon, "Closed Eyes"
        'E': ('French', 1906),  # Georges Rouault, "The Drinker"
        'F': ('N/A', 2023)      # AI-generated image
    }
    
    result_string = ""
    
    # Iterate through the paintings in alphabetical order of their IDs
    for painting_id in sorted(paintings.keys()):
        nationality, year = paintings[painting_id]
        
        # Check if the painter is French and the year is before 1900
        if nationality == 'French' and year < 1900:
            result_string += 'Y'
        else:
            result_string += 'N'
            
    print(result_string)

solve_art_puzzle()