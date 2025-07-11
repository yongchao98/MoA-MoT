def solve_art_puzzle():
    """
    This function analyzes six artworks to determine which were created by a French painter before 1900.
    """
    
    # Information about each artwork:
    # A: Félix Vallotton (Swiss until 1900), 'The Bathing on a Summer Evening', 1892-93. Condition: Not French at time of creation.
    # B: Édouard Vuillard (French), 'Le dépôt. Chez le brocanteur', c. 1892. Condition: Meets all criteria.
    # C: Henri-Edmond Cross (French), 'La Ferme, Matin', 1893. Condition: Meets all criteria.
    # D: Odilon Redon (French), 'Closed Eyes', 1890. Condition: Meets all criteria.
    # E: Henri de Toulouse-Lautrec (French), 'The Hangover', c. 1888. Condition: Meets all criteria.
    # F: Modern digital/AI creation. Condition: Not from before 1900.
    
    artworks = [
        {'id': 'A', 'is_french': False, 'before_1900': True},
        {'id': 'B', 'is_french': True,  'before_1900': True},
        {'id': 'C', 'is_french': True,  'before_1900': True},
        {'id': 'D', 'is_french': True,  'before_1900': True},
        {'id': 'E', 'is_french': True,  'before_1900': True},
        {'id': 'F', 'is_french': False, 'before_1900': False} # Considered modern, so fails both
    ]
    
    result_string = ""
    for art in artworks:
        if art['is_french'] and art['before_1900']:
            result_string += "Y"
        else:
            result_string += "N"
            
    print(result_string)

solve_art_puzzle()