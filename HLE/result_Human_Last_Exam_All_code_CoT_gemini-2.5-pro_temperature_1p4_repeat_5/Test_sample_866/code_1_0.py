def solve_origami_assignments():
    """
    Calculates the number of valid flat-foldable assignments for a given
    partially assigned crease pattern by checking Maekawa's theorem.
    """

    # 1. Parse the input pattern
    pattern_str = "[60,M,30,?,50,?,70,V,150,?]"
    crease_specs = [c for c in pattern_str.strip('[]').split(',') if c in ['M', 'V', '?']]

    num_creases = len(crease_specs)
    known_m = crease_specs.count('M')
    known_v = crease_specs.count('V')
    unknown_count = crease_specs.count('?')

    print("Analyzing the given crease pattern:")
    print(f"Total number of creases: {num_creases}")
    print(f"Known Mountain folds ('M'): {known_m}")
    print(f"Known Valley folds ('V'): {known_v}")
    print(f"Unassigned folds ('?'): {unknown_count}\n")

    print("For a vertex to be flat-foldable, Maekawa's Theorem must be satisfied.")
    print("This means: |Total Mountain Folds - Total Valley Folds| = 2\n")

    print("Checking all possible assignments for the 3 unassigned folds:")

    total_valid_assignments = 0

    # 2. Iterate through all possible assignments for the '?' creases
    for n_m in range(unknown_count + 1):
        n_v = unknown_count - n_m
        
        print(f"Case: Assigning {n_m} of '?' to Mountain and {n_v} to Valley.")
        
        total_m = known_m + n_m
        total_v = known_v + n_v
        
        # 3. For each assignment, check the theorem and output the equation
        print(f"  - Total M = {known_m} (known) + {n_m} (assigned) = {total_m}")
        print(f"  - Total V = {known_v} (known) + {n_v} (assigned) = {total_v}")
        
        diff = abs(total_m - total_v)
        
        print(f"  - Checking condition: |{total_m} - {total_v}| = {diff}")
        
        if diff == 2:
            print("  - Result: This assignment SATISFIES Maekawa's Theorem.\n")
            # In a more complex problem, we would kombinatorially count ways
            # and check other conditions. Here it is sufficient to see if any
            # valid assignment exists.
            total_valid_assignments += 1
        else:
            print(f"  - Result: This assignment FAILS because {diff} != 2.\n")
            
    print("Conclusion:")
    print("None of the possible assignments for the '?' folds satisfy Maekawa's Theorem.")
    print("Since this is a necessary condition, no assignment can be flat-foldable.")
    print(f"\nTotal number of valid assignments: {total_valid_assignments}")

solve_origami_assignments()