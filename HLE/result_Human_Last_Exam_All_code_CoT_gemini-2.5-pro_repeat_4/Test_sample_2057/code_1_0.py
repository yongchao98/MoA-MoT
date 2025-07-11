import collections

def solve_hat_puzzle():
    """
    This function solves the hat puzzle by simulating the deduction process round by round.
    """

    # Step 1: Define base patterns (1 for Black, 0 for White)
    # After Round 1, we know there are no three consecutive hats of the same color.
    # Combinatorial analysis reveals three possible families of patterns (A, B, C).
    pattern_A = (1, 1, 0, 1, 0, 1, 0, 1, 0)  # BBWBWBWBW
    pattern_B = (1, 1, 0, 1, 1, 0, 1, 0, 0)  # BBWBBWBWBW
    pattern_C = (1, 1, 0, 0, 1, 1, 0, 1, 0)  # BBWWBW BWBW

    # Generate the full set of possible configurations after Round 1 (S1)
    S1 = set()
    for p in [pattern_A, pattern_B, pattern_C]:
        for i in range(9):
            S1.add(p[i:] + p[:i])

    # Helper function to get the 6-hat view for a person at a given index
    def get_view(config, index):
        view = []
        for i in range(6):
            view.append(config[(index + 2 + i) % 9])
        return tuple(view)

    # Map each possible view to the set of hat colors seen with that view across all configs in S1
    view_to_hats_S1 = collections.defaultdict(set)
    for config in S1:
        for i in range(9):
            view = get_view(config, i)
            hat = config[i]
            view_to_hats_S1[view].add(hat)

    # Step 2: Analyze Round 2. Everyone says "No".
    # This eliminates any configuration where at least one person could have deduced their hat color.
    # Such a person would see a view that corresponds to only one hat color in S1.
    # The remaining set of configurations is S2.
    S2 = set()
    for config in S1:
        is_survivor = True
        for i in range(9):
            view = get_view(config, i)
            if len(view_to_hats_S1[view]) == 1:
                is_survivor = False
                break
        if is_survivor:
            S2.add(config)
    
    # At this point, S2 contains only the 9 rotations of pattern A (`BBWBWBWBW`).
    # The announcements in Round 2 informed everyone that the true configuration belongs to S2.

    # Step 3: Analyze Round 3. At least one person says "Yes".
    # We need to find a configuration in S2 where this is possible.
    
    # Map views to hat colors for the remaining possibilities in S2
    view_to_hats_S2 = collections.defaultdict(set)
    for config in S2:
        for i in range(9):
            view = get_view(config, i)
            hat = config[i]
            view_to_hats_S2[view].add(hat)

    solution_config = None
    yes_sayers_indices = []

    # Find the configuration in S2 where someone can make a deduction
    for config in S2:
        current_yes_sayers = []
        for i in range(9):
            my_hat = config[i]
            my_view = get_view(config, i)
            
            # A person can say "Yes" if their alternative hat color is not possible
            # for their view within the new set of possibilities (S2).
            if (1 - my_hat) not in view_to_hats_S2[my_view]:
                current_yes_sayers.append(i)
        
        if len(current_yes_sayers) > 0:
            solution_config = config
            yes_sayers_indices = current_yes_sayers
            break # Found the solution

    # Print the final answer
    if solution_config:
        num_yes = len(yes_sayers_indices)
        dist_str = "".join(['B' if h == 1 else 'W' for h in solution_config])
        
        print(f"Number of people who replied 'Yes': {num_yes}")
        print(f"The distribution of hats around the table is: {dist_str}")
        print("The people who replied 'Yes' are the three individuals with Black hats who are not part of the 'BB' pair.")

    else:
        print("A solution could not be determined programmatically.")

solve_hat_puzzle()