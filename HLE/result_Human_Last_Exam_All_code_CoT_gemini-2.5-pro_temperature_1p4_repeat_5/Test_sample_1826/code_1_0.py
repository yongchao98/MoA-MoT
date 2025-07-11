def solve():
    """
    This function demonstrates the construction of an infinite set x that is
    almost disjoint from every set in a countable family S of infinite sets.
    """
    # Let S be a countable family of infinite subsets of omega.
    # For demonstration, we use a few example sets.
    s0 = {n for n in range(0, 100, 2)}   # Even numbers
    s1 = {n for n in range(0, 100, 3)}   # Multiples of 3
    s2 = {n for n in range(0, 100) if pow(n, 0.5).is_integer()} # Perfect squares
    s3 = {2**n -1 for n in range(1, 7)} # Mersenne numbers
    
    S = [s0, s1, s2, s3]
    for i in range(len(S)):
        # Sort the sets to make finding the minimum element easier
        S[i] = sorted(list(S[i]))
    
    x = []
    x_prev = -1
    num_elements_to_construct = 20

    print("Constructing the set x...")
    
    for n in range(num_elements_to_construct):
        # At step n, we consider sets s_i for i < n.
        
        # The set of values to "jump over".
        # It contains x_{n-1} and the smallest element > x_{n-1} from each s_i (for i<n).
        values_to_surpass = [x_prev]
        
        # Consider the first n sets from S (or fewer if S is smaller than n).
        num_sets_to_consider = min(n, len(S))

        for i in range(num_sets_to_consider):
            # Find the smallest element in s_i that is greater than x_prev
            y_i = -1 
            for val in S[i]:
                if val > x_prev:
                    y_i = val
                    break
            
            # If s_i has no element greater than x_prev in our sample, we can ignore it.
            # In theory, s_i is infinite, so such an element always exists.
            if y_i != -1:
                values_to_surpass.append(y_i)
        
        # x_n is one greater than the maximum of these values.
        x_n = max(values_to_surpass) + 1
        
        print(f"Step {n}:")
        print(f"  x_{n-1} = {x_prev}")
        if n > 0:
            print(f"  Values to surpass: {values_to_surpass}")
        print(f"  Chose x_{n} = {x_n}")

        x.append(x_n)
        x_prev = x_n

    print("\nConstructed set x = ", sorted(x))
    
    # Verify the property for the constructed part of x
    print("\nVerifying intersections:")
    for i in range(len(S)):
        intersection = set(x).intersection(set(S[i]))
        print(f"|x intersect s{i}| = {len(intersection)}. Intersection: {sorted(list(intersection))}")

solve()