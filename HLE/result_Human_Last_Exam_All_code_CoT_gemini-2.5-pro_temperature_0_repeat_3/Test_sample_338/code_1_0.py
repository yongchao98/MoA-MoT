import itertools

def get_hat_type(hat):
    """Returns 'C' for color hats and 'N' for number hats."""
    if isinstance(hat, str):
        return 'C'
    return 'N'

def knows(person_hat, seen_hats, all_hats):
    """
    Determines if a person can deduce their hat based on what they see.
    """
    n_b = all_hats.count('B')
    n_w = all_hats.count('W')
    n_n = 7 - n_b - n_w

    seen_b = seen_hats.count('B')
    seen_w = seen_hats.count('W')
    seen_n = sum(1 for h in seen_hats if isinstance(h, int))

    # Case 1: Deduces color by absence.
    if n_b > 0 and seen_b == 0:
        return True
    if n_w > 0 and seen_w == 0:
        return True

    # Case 2: Deduces it's the only number hat.
    if n_n == 1 and seen_n == 0:
        return True

    # Case 3: Deduces type is Number and knows the exact number from a gap.
    # A person knows their type is Number if they see at least one of each color.
    if seen_b > 0 and seen_w > 0:
        if isinstance(person_hat, int):
            seen_numbers = sorted([h for h in seen_hats if isinstance(h, int)])
            
            # Not enough numbers to find a gap.
            if len(seen_numbers) < 1:
                return False

            # Check for a single gap of size 2 (e.g., sees 3, 5 -> must be 4)
            # This is the only way to uniquely determine the number from a gap.
            gaps = []
            for i in range(len(seen_numbers) - 1):
                diff = seen_numbers[i+1] - seen_numbers[i]
                if diff > 1:
                    gaps.append(diff)
            
            if len(gaps) == 1 and gaps[0] == 2:
                return True
    
    return False

def check_config(config):
    """
    Checks if a given hat configuration matches the K/DK sequence.
    """
    expected_results = [True, False, True, False, True, False, True] # K, DK, K, DK, K, DK, K
    
    for i in range(7):
        person_hat = config[i]
        seen_hats = config[:i] + config[i+1:]
        result = knows(person_hat, seen_hats, config)
        if result != expected_results[i]:
            return False
    return True

def solve_puzzle():
    """
    Finds all valid configurations and identifies Alice.
    """
    people = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    valid_type_configs = set()

    # Iterate through all possible hat counts (n_b, n_w, n_n)
    for n_b in range(1, 6):
        for n_w in range(1, 6):
            n_n = 7 - n_b - n_w
            if n_n < 1:
                continue

            base_hats = ['B'] * n_b + ['W'] * n_w + list(range(1, n_n + 1))
            
            # Iterate through all unique permutations of these hats
            for p in set(itertools.permutations(base_hats)):
                config = list(p)
                if check_config(config):
                    type_config = tuple(get_hat_type(h) for h in config)
                    valid_type_configs.add(type_config)

    # Find Alice by comparing the valid type configurations
    valid_list = list(valid_type_configs)
    ambiguous_indices = set()

    if len(valid_list) < 2:
        print("Could not find two distinct valid scenarios.")
        return

    for i in range(len(valid_list)):
        for j in range(i + 1, len(valid_list)):
            t1 = valid_list[i]
            t2 = valid_list[j]
            
            diff_indices = [k for k in range(7) if t1[k] != t2[k]]
            
            if len(diff_indices) == 1:
                ambiguous_indices.add(diff_indices[0])

    if len(ambiguous_indices) == 1:
        alice_index = list(ambiguous_indices)[0]
        alice_name = people[alice_index]
        print(f"Found two scenarios that produce the same result.")
        print(f"These scenarios differ only on the hat type of person {alice_name}.")
        print(f"Therefore, Alice is {alice_name}.")
        print(f"<<<{alice_name}>>>")
    else:
        print("Could not uniquely identify Alice based on the logic.")

solve_puzzle()