import collections
from itertools import combinations

def find_ca_rules():
    """
    Analyzes the evolution of a 1D cellular automaton to find all possible
    elementary rules that could have generated it.
    """
    # 1. Represent the pattern from the image as a grid of 0s (white) and 1s (black).
    grid = [
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0],
        [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
    ]
    height = len(grid)
    width = len(grid[0])

    # 2. Extract all observed transitions (neighborhood -> next state).
    observed_transitions = {}
    for t in range(height - 1):
        for i in range(width):
            # Assume cells outside the grid are 0 (fixed boundary condition).
            left = grid[t][i - 1] if i > 0 else 0
            center = grid[t][i]
            right = grid[t][i + 1] if i < width - 1 else 0
            
            neighborhood = (left, center, right)
            next_state = grid[t + 1][i]

            # Store the observed transition if it's new.
            # A valid ECA will not have conflicting rules.
            if neighborhood not in observed_transitions:
                observed_transitions[neighborhood] = next_state

    # The standard Wolfram ordering for elementary CA rules.
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    # Map each neighborhood to its corresponding bit position in the rule number.
    neighborhood_to_bit_pos = {n: 7 - i for i, n in enumerate(neighborhood_order)}

    # 3. Identify unobserved neighborhoods and calculate a base rule.
    base_rule = 0
    wildcard_powers = []
    
    for neighborhood in neighborhood_order:
        bit_pos = neighborhood_to_bit_pos[neighborhood]
        if neighborhood in observed_transitions:
            # If the transition is observed, its output bit is fixed.
            bit_value = observed_transitions[neighborhood]
            if bit_value == 1:
                base_rule += 2**bit_pos
        else:
            # If unobserved, its output can be 0 or 1. This is a "wildcard".
            # We store the value it would add to the rule if the output were 1.
            wildcard_powers.append(2**bit_pos)

    # 4. Generate all possible rules by adding combinations of wildcard values.
    possible_rules = set()
    num_wildcards = len(wildcard_powers)

    for i in range(num_wildcards + 1):
        for combo in combinations(wildcard_powers, i):
            rule = base_rule + sum(combo)
            possible_rules.add(rule)
    
    # 5. Sort the resulting rules and print them in the required format.
    sorted_rules = sorted(list(possible_rules))
    print(','.join(map(str, sorted_rules)))

find_ca_rules()