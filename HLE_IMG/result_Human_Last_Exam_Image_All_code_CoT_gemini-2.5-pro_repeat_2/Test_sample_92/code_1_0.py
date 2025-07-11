def solve_art_puzzle():
    """
    This function determines which paintings were created by a French painter before 1900.
    """
    paintings_info = [
        {'id': 'A', 'artist': 'Félix Vallotton', 'nationality': 'Swiss', 'year': 1897},
        {'id': 'B', 'artist': 'Édouard Manet', 'nationality': 'French', 'year': 1881},
        {'id': 'C', 'artist': 'Henri-Edmond Cross', 'nationality': 'French', 'year': 1895},
        {'id': 'D', 'artist': 'Odilon Redon', 'nationality': 'French', 'year': 1890},
        {'id': 'E', 'artist': 'Georges Rouault', 'nationality': 'French', 'year': 1907},
        {'id': 'F', 'artist': 'Unknown (Modern/AI)', 'nationality': 'N/A', 'year': 2000} # Placeholder for modern art
    ]

    result_string = ""
    for painting in paintings_info:
        # Check if the painter was French and the work was created before 1900
        is_french_before_1900 = painting['nationality'] == 'French' and painting['year'] < 1900
        
        if is_french_before_1900:
            result_string += 'Y'
        else:
            result_string += 'N'
            
    print("The final answer is:")
    print(result_string)

solve_art_puzzle()