import collections
import itertools

def get_view(config, person_index):
    """Gets the 6 hats a person can see."""
    n = len(config)
    view = []
    for i in range(3, n):
        view.append(config[(person_index + i) % n])
    return tuple(sorted(view))

def could_deduce_r1(view):
    """Checks if a view allows a deduction in Round 1."""
    # 5 Black hats seen -> own must be White
    # 4 White hats seen -> own must be Black
    counts = collections.Counter(view)
    if counts['B'] == 5 or counts['W'] == 4:
        return True
    return False

def find_solution():
    """
    Solves the hat puzzle by simulating the rounds of deductions.
    """
    hats = ['B'] * 5 + ['W'] * 4
    # Generate unique configurations (handling rotations)
    all_perms = set(itertools.permutations(hats))
    
    unique_configs = []
    seen_configs = set()
    for p in all_perms:
        is_new = True
        temp_p = collections.deque(p)
        for _ in range(len(p)):
            if tuple(temp_p) in seen_configs:
                is_new = False
                break
            temp_p.rotate(1)
        if is_new:
            unique_configs.append(p)
            # Add all rotations to seen to avoid duplicates
            temp_p_2 = collections.deque(p)
            for _ in range(len(p)):
                seen_configs.add(tuple(temp_p_2))
                temp_p_2.rotate(1)

    # --- Round 1 ---
    possible_after_r1 = []
    for config in unique_configs:
        someone_deduced = False
        for i in range(len(config)):
            view = get_view(config, i)
            if could_deduce_r1(view):
                someone_deduced = True
                break
        if not someone_deduced:
            possible_after_r1.append(config)

    # --- Round 2 ---
    possible_after_r2 = []
    for config in possible_after_r1:
        someone_deduced = False
        for i in range(len(config)):
            my_view = get_view(config, i)
            my_hat = config[i]
            
            # Hypothesis: my hat is the OTHER color
            other_color = 'W' if my_hat == 'B' else 'B'
            
            # Is there any world in possible_after_r1 consistent with my view and the other hat color?
            is_other_possible = False
            for other_config in possible_after_r1:
                if get_view(other_config, i) == my_view and other_config[i] == other_color:
                    is_other_possible = True
                    break
            
            # If the other color is not possible, I can deduce my color.
            if not is_other_possible:
                someone_deduced = True
                break
        if not someone_deduced:
            possible_after_r2.append(config)

    # --- Round 3 ---
    final_configs = []
    num_yes = 0
    for config in possible_after_r2:
        can_deduce_in_r3 = False
        deducers = []
        for i in range(len(config)):
            my_view = get_view(config, i)
            my_hat = config[i]
            other_color = 'W' if my_hat == 'B' else 'B'
            
            is_other_possible = False
            for other_config in possible_after_r2:
                if get_view(other_config, i) == my_view and other_config[i] == other_color:
                    is_other_possible = True
                    break
            
            if not is_other_possible:
                can_deduce_in_r3 = True
                deducers.append(i)

        if can_deduce_in_r3:
            final_configs.append(config)
            num_yes = len(deducers)

    if final_configs:
        # There should be only one unique configuration family as the answer
        solution_config = final_configs[0]
        print(f"Number of people who replied 'Yes': {num_yes}")
        print("The distribution of hats around the table is:")
        print(' '.join(solution_config))
    else:
        print("No solution found.")

find_solution()