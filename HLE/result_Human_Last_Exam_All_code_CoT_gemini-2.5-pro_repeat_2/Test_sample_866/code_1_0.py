import itertools

def solve_crease_pattern():
    """
    Calculates the number of valid flat-fordable assignments for a given crease pattern.
    """
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
    
    angles = [p for p in pattern if isinstance(p, int)]
    creases = [p for p in pattern if isinstance(p, str)]

    # Step 1: Check Kawasaki's Theorem
    # The sum of angles around a vertex must be 360 degrees.
    # They must also be partitionable into two sets that each sum to 180.
    if sum(angles) != 360:
        print("Sum of angles is not 360. Not flat-fordable.")
        print(0)
        return

    # Find a subset of angles that sums to 180
    found_partition = False
    n = len(angles)
    for i in range(1, 1 << n):
        subset1 = []
        subset2 = []
        current_sum = 0
        for j in range(n):
            if (i >> j) & 1:
                current_sum += angles[j]
                subset1.append(angles[j])
        
        if current_sum == 180:
            # Found a valid partition.
            for j in range(n):
                if not ((i >> j) & 1):
                    subset2.append(angles[j])
            
            print("Kawasaki's Theorem is satisfied. The angle partitions are:")
            print(f"{' + '.join(map(str, subset1))} = 180")
            print(f"{' + '.join(map(str, subset2))} = 180")
            found_partition = True
            break
            
    if not found_partition:
        print("Kawasaki's Theorem is not satisfied.")
        print(0)
        return

    # Step 2: Check Maekawa's Theorem for all possible assignments
    num_unknowns = creases.count('?')
    unknown_indices = [i for i, c in enumerate(creases) if c == '?']
    
    valid_assignments_count = 0
    
    # Iterate through all 2^num_unknowns possibilities for '?'
    for p in itertools.product(['M', 'V'], repeat=num_unknowns):
        temp_creases = list(creases)
        for i, val in enumerate(p):
            temp_creases[unknown_indices[i]] = val
        
        num_m = temp_creases.count('M')
        num_v = temp_creases.count('V')
        
        # Check Maekawa's Theorem: |#M - #V| = 2
        if abs(num_m - num_v) == 2:
            valid_assignments_count += 1
            # If any valid assignment was found, we would print its equation here.
            # e.g., print(f"Valid assignment found: {temp_creases}. Equation: {num_m} - {num_v} = {num_m - num_v}")

    if valid_assignments_count == 0:
        print("\nNo assignment satisfies Maekawa's Theorem (|#M - #V| = 2).")

    print(f"\nTotal number of different flat-fordable assignments: {valid_assignments_count}")

solve_crease_pattern()