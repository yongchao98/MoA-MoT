import collections
from itertools import permutations

def solve_hat_puzzle():
    """
    Solves the 7-person hat puzzle to identify Alice.
    
    This function generates all possible valid hat configurations, simulates
    the dialogue to filter these configurations, and then analyzes the
    surviving scenarios to determine for which person the hat type
    (Color vs. Number) remains ambiguous.
    """
    
    # Hat representation: Black=-1, White=-2, Numbers>=0
    
    def get_hat_type(hat_val):
        return "Color" if hat_val < 0 else "Number"

    # Step 1: Generate all possible valid configurations
    all_configs = set()
    # n_N: number of Number hats
    # n_B: number of Black hats
    # n_W: number of White hats
    for n_N in range(1, 6):  # 1 to 5 number hats
        numbers = list(range(n_N))
        for n_B in range(1, 7 - n_N):
            n_W = 7 - n_N - n_B
            if n_W < 1:
                continue

            hat_pool = [-1] * n_B + [-2] * n_W + numbers
            
            # Use set to get unique permutations
            for p in set(permutations(hat_pool)):
                all_configs.add(p)

    initial_possible_worlds = tuple(all_configs)

    # Step 2: Simulate the dialogue to filter the configurations
    dialogue = [True, False, True, False, True, False, True]  # True="I know", False="I don't know"
    person_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    
    current_possible_worlds = initial_possible_worlds

    for i in range(7):  # For each person A through G
        person_knows_statement = dialogue[i]
        
        # Group worlds by what the current person (i) would see
        views = collections.defaultdict(list)
        for world in current_possible_worlds:
            # The view is the configuration tuple without the person's own hat
            view = world[:i] + world[i+1:]
            views[view].append(world)

        next_possible_worlds = []
        for view, worlds_with_this_view in views.items():
            
            # Determine if the person would know their hat given this view.
            # They know if only one hat is possible for them, OR if they see
            # no number hats (which forces their own hat to be a number).
            possible_hats = {world[i] for world in worlds_with_this_view}
            view_has_no_numbers = all(hat < 0 for hat in view)

            would_know = (len(possible_hats) == 1) or view_has_no_numbers
            
            # If the person's deduction matches the dialogue, these worlds survive
            if would_know == person_knows_statement:
                next_possible_worlds.extend(worlds_with_this_view)
        
        current_possible_worlds = tuple(next_possible_worlds)

    # Step 3: Analyze the surviving configurations to find Alice
    surviving_configs = current_possible_worlds
    
    if not surviving_configs:
        print("No solution found. The puzzle parameters may lead to a contradiction.")
        return

    alice_candidates = []
    for i in range(7):
        # Check if the hat type for person i is consistent across all solutions
        hat_types_for_person_i = {get_hat_type(world[i]) for world in surviving_configs}
        
        # If the type is not consistent (i.e., sometimes Color, sometimes Number), this person is Alice.
        if len(hat_types_for_person_i) > 1:
            alice_candidates.append(person_names[i])

    # Final Output
    if len(alice_candidates) == 1:
        alice = alice_candidates[0]
        print(f"The blind person H can determine the hat type for everyone except for one person.")
        print(f"The person whose hat type (Color or Number) is ambiguous is: {alice}")
        print(f"Therefore, Alice is {alice}.")
        print(f"<<<{alice}>>>")
    elif len(alice_candidates) == 0:
        print("A unique solution for Alice was not found: H can determine everyone's hat type.")
    else:
        print(f"A unique solution for Alice was not found: H cannot determine hat types for multiple people: {alice_candidates}")


solve_hat_puzzle()