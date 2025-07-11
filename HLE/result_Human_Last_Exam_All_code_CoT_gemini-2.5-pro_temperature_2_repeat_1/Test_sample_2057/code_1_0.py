import collections
from itertools import combinations, permutations

def solve_hat_puzzle():
    """
    Solves the 9-person hat puzzle by iteratively filtering possible arrangements.
    """
    # --- Configuration ---
    NUM_PEOPLE = 9
    TOTAL_HATS = {'B': 5, 'W': 4}

    # --- Helper Functions ---
    def get_canonical(s):
        """Gets the lexicographically smallest rotation of a circular arrangement."""
        return min(s[i:] + s[:i] for i in range(len(s)))

    def generate_all_arrangements():
        """Generates all C(9,4) unique hat arrangements."""
        base = ['B'] * NUM_PEOPLE
        arrangements = []
        for indices in combinations(range(NUM_PEOPLE), TOTAL_HATS['W']):
            arr = list(base)
            for i in indices:
                arr[i] = 'W'
            arrangements.append("".join(arr))
        return arrangements

    def get_seen_hats(arrangement, person_index):
        """Gets the 6 hats a person at a given index sees."""
        n = len(arrangement)
        seen_hats_list = []
        for i in range(2, n - 1):
            idx = (person_index + i) % n
            seen_hats_list.append(arrangement[idx])
        return "".join(seen_hats_list)

    # --- Deduction Logic Functions ---
    def can_deduce_r1(seen_hats):
        """Checks for the simple Round 1 deduction."""
        if seen_hats.count('W') == TOTAL_HATS['W']:
            return 'B'  # All white hats seen, so unseen must be black.
        if seen_hats.count('B') == TOTAL_HATS['B']:
            return 'W'  # All black hats seen, so unseen must be white.
        return None

    def is_r1_valid(arrangement):
        """Checks if an arrangement survives Round 1 (no one says 'Yes')."""
        for i in range(NUM_PEOPLE):
            if can_deduce_r1(get_seen_hats(arrangement, i)):
                return False
        return True

    def can_deduce(seen_hats, person_index, possible_worlds_set):
        """
        Generic deduction logic for Rounds 2 and 3.
        A person can deduce if their hat color is the same across all possible
        worlds that are consistent with what they see.
        """
        # Filter the possible worlds to those consistent with the person's view
        consistent_worlds = [
            world for world in possible_worlds_set 
            if get_seen_hats(world, person_index) == seen_hats
        ]
        
        if not consistent_worlds:
            return None # Should not happen in this puzzle's logic flow.

        # Check the person's hat color in all these consistent worlds
        possible_hat_colors = {world[person_index] for world in consistent_worlds}

        if len(possible_hat_colors) == 1:
            return list(possible_hat_colors)[0] # Deduction is successful
        
        return None # Still ambiguous

    # --- Main Solving Logic ---

    # Step 1: Generate all arrangements
    all_arrangements = generate_all_arrangements()

    # Step 2: Filter for Round 1 survivors
    r1_valid_set = set()
    for arr in all_arrangements:
        if is_r1_valid(arr):
            # Add all rotations of this valid arrangement to our set
            for i in range(NUM_PEOPLE):
                r1_valid_set.add(arr[i:] + arr[:i])

    # Step 3: Filter for Round 2 survivors
    r2_valid_set = set()
    for arr in r1_valid_set:
        is_deducible_in_r2 = False
        for p_idx in range(NUM_PEOPLE):
            seen_by_p = get_seen_hats(arr, p_idx)
            if can_deduce(seen_by_p, p_idx, r1_valid_set):
                is_deducible_in_r2 = True
                break
        if not is_deducible_in_r2:
            r2_valid_set.add(arr)

    # Step 4: Analyze Round 3 to find the final answer
    solution = None
    for arr in r2_valid_set:
        num_yes_in_r3 = 0
        for p_idx in range(NUM_PEOPLE):
            seen_by_p = get_seen_hats(arr, p_idx)
            if can_deduce(seen_by_p, p_idx, r2_valid_set):
                num_yes_in_r3 += 1
        
        if num_yes_in_r3 > 0:
            # Found the solution that matches the puzzle description
            solution = {"arrangement": arr, "count": num_yes_in_r3}
            break

    # Step 5: Print the results
    if solution:
        num_people = solution['count']
        arrangement_str = " ".join(list(solution['arrangement']))
        
        print(f"The number of people who replied 'Yes' is: {num_people}")
        print(f"The distribution of hats around the table is: {arrangement_str}")
        print("(B = Black, W = White; other rotations are also possible)")
    else:
        print("No solution was found based on the puzzle's logic.")

solve_hat_puzzle()
<<<4>>>