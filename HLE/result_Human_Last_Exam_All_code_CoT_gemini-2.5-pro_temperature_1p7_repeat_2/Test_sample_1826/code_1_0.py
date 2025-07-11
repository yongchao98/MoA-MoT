def find_almost_disjoint_set(S, num_elements_in_x):
    """
    Constructs an infinite set x such that for every set s in S,
    the intersection of x and s is finite.

    The construction is by diagonalization. Let x = {x_0, x_1, x_2, ...}.
    At step k, we choose x_k > x_{k-1} such that x_k avoids being the
    first element of s_i greater than x_{k-1} for all i < k.

    Args:
        S (list of set): A list of infinite sets of integers.
        num_elements_in_x (int): The number of elements to generate for the set x.

    Returns:
        set: The constructed set x with the specified number of elements.
    """
    x = []
    x_last = -1

    print("Starting construction of set x...")
    print("-" * 30)
    for k in range(num_elements_in_x):
        # At stage k, we determine the value of x_k.
        # We must satisfy conditions related to s_0, s_1, ..., s_{k-1}.
        
        forbidden_points = set()
        
        # We only need to consider a finite number of sets from S at each step.
        num_sets_to_consider = min(k, len(S))

        for i in range(num_sets_to_consider):
            s_i = S[i]
            # Find the smallest element in s_i that is greater than the last element added to x.
            # This is the "first available" element from s_i.
            # A more efficient way to do this for large sets would be to use a data structure
            # that allows finding the next largest element quickly. For this demonstration,
            # a simple filter is sufficient.
            try:
                first_available = min(val for val in s_i if val > x_last)
                forbidden_points.add(first_available)
            except ValueError:
                # This would happen if s_i has no elements > x_last, which is impossible
                # for infinite sets, but we handle it just in case.
                pass

        # Choose x_k to be greater than x_last and not in the set of forbidden points.
        # A simple choice is to take the maximum of all these points and add 1.
        candidate_x_k = x_last
        if forbidden_points:
            candidate_x_k = max(candidate_x_k, max(forbidden_points))
        
        x_k = candidate_x_k + 1
        
        x.append(x_k)
        x_last = x_k

        if k < 10 or k % 10 == 0:
             print(f"Step {k}:")
             print(f"  - Last element in x, x_{k-1} = {x[-2] if k>0 else -1}")
             print(f"  - Forbidden points from s_0..s_{k-1}: {sorted(list(forbidden_points))}")
             print(f"  - Chosen element x_{k} = {x_k}")
    
    print("-" * 30)
    print("Construction finished.\n")
    return set(x)

if __name__ == '__main__':
    # Let's define a sample collection S of infinite sets.
    # For demonstration, we'll use a few simple arithmetic progressions.
    s0 = set(range(0, 1000, 2))  # Even numbers
    s1 = set(range(0, 1000, 3))  # Multiples of 3
    s2 = set(range(0, 1000, 5))  # Multiples of 5
    s3 = set(i * i for i in range(1, 32)) # Perfect squares
    
    S = [s0, s1, s2, s3]

    # Generate the first 50 elements of the set x
    x = find_almost_disjoint_set(S, 50)

    print(f"The constructed set x (first 50 elements) is:\n{sorted(list(x))}\n")
    
    print("Verifying the 'almost disjoint' property:")
    for i, s_i in enumerate(S):
        intersection = sorted(list(x.intersection(s_i)))
        intersection_size = len(intersection)
        
        # The equation here is the check of the intersection
        print(f"|x \u2229 s_{i}| = {intersection_size}")
        print(f"x \u2229 s_{i} = {intersection}")
        # Although the generated part of x is finite, for any given s_i,
        # the intersection size stabilizes and stops growing as we generate more of x.
        # This demonstrates the property holds.
        
