def construct_almost_disjoint_set(S, num_elements):
    """
    Constructs the first `num_elements` of an infinite set x that is
    almost disjoint from every set in the countable collection S.

    Args:
        S (list of list of int): A list of lists, where each inner list
                                 represents an infinite set of natural numbers,
                                 sorted in increasing order.
        num_elements (int): The number of elements to construct for the set x.

    Returns:
        list of int: The first `num_elements` of the constructed set x.
    """
    x = []
    x_prev = -1

    print("--- Construction of set x ---")

    for n in range(num_elements):
        # We need to find max(x_{n-1}, a_{0,n}, a_{1,n}, ..., a_{n,n})
        # where a_{i,n} is the n-th element of set s_i.
        
        # Ensure we don't go out of bounds for the number of sets in S
        # or the number of elements in the example sets.
        max_n = min(n, len(S) - 1)
        
        candidate_values = [x_prev]
        print(f"Step n={n}:")
        print(f"  Previous element x_{n-1} = {x_prev}")
        
        for i in range(max_n + 1):
            s_i = S[i]
            # Make sure the example sets are long enough for the demonstration
            if n < len(s_i):
                a_in = s_i[n]
                candidate_values.append(a_in)
                print(f"  From s_{i}, taking element a_{i},{n} = {a_in}")
            else:
                print(f"  s_{i} is not long enough to get element at index {n}. Skipping.")

        # Calculate the new element for x
        m_n = max(candidate_values)
        x_n = m_n + 1
        x.append(x_n)
        
        print(f"  Candidate values for max: {candidate_values}")
        print(f"  max = {m_n}")
        print(f"  Constructed x_{n} = {m_n} + 1 = {x_n}\n")
        
        x_prev = x_n
        
    return x

def main():
    # Define a sample countable collection of infinite sets S.
    # s_0: even numbers
    # s_1: multiples of 3
    # s_2: multiples of 5
    # s_3: powers of 2 (starting from 2^1)
    s0 = [i * 2 for i in range(100)]
    s1 = [i * 3 for i in range(100)]
    s2 = [i * 5 for i in range(100)]
    s3 = [2**i for i in range(1, 101)]
    S = [s0, s1, s2, s3]

    num_x_elements_to_construct = 10
    x = construct_almost_disjoint_set(S, num_x_elements_to_construct)

    print("--- Verification ---")
    print(f"Constructed set x = {x}")

    for i, s_i in enumerate(S):
        # Convert to sets for efficient intersection
        intersection = set(x).intersection(set(s_i))
        print(f"Intersection of x with s_{i}: |x intersect s_{i}| = |{list(intersection)}| = {len(intersection)}")
        # The proof guarantees the intersection is finite. Our small sample should reflect this.
        
if __name__ == "__main__":
    main()