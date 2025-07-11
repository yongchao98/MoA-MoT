import itertools

def find_largest_matchup_list():
    """
    This function solves the Tennis matchup problem by using a greedy algorithm.
    It finds a large set of 4-player matchups from 11 players such that
    no two matchups have more than 2 players in common.
    """
    
    # Step 1: Define the players. We use numbers 0 to 10 for simplicity.
    players = list(range(11))
    
    # Step 2: Generate all possible matchups (4-player combinations).
    all_matchups = list(itertools.combinations(players, 4))
    
    # Step 3: Use a greedy approach to build the list of valid matchups.
    selected_matchups = []
    
    # Iterate through all generated matchups in their default lexicographical order.
    for candidate_matchup in all_matchups:
        candidate_set = set(candidate_matchup)
        is_valid = True
        
        # Check the candidate against all previously selected matchups.
        for existing_matchup in selected_matchups:
            existing_set = set(existing_matchup)
            
            # Calculate how many players are common to both matchups.
            intersection_size = len(candidate_set.intersection(existing_set))
            
            # If the intersection is greater than 2, the candidate is invalid.
            if intersection_size > 2:
                is_valid = False
                break
        
        # If the candidate is valid, add it to our list of selected matchups.
        if is_valid:
            selected_matchups.append(candidate_matchup)
            
    # Step 4: Output the final result in the requested format.
    num_matchups = len(selected_matchups)
    print(f"The largest list of matchups found has a size of: {num_matchups}")
    print("The final equation for the number of matchups is:")
    
    # Create the equation string, e.g., "1 + 1 + ... + 1 = total"
    equation_parts = ["1"] * num_matchups
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {num_matchups}")

# Execute the function to find and print the solution.
find_largest_matchup_list()
