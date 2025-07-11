import math

def cantor_pairing(k1, k2):
    """
    Computes the Cantor pairing function for two non-negative integers.
    """
    return (k1 + k2) * (k1 + k2 + 1) // 2 + k2

def construct_f_S(S):
    """
    Returns a function f_S(n) = pi(n, |S intersect {1, ..., n}|).
    """
    # For efficiency in the demo, convert S to a set data structure.
    S_set = set(S)
    memoized_intersection_size = {0: 0}

    def f_S(n):
        """
        Computes the value of the function f_S at n.
        This implementation is for demonstration and not optimized for very large n.
        """
        # Calculate |S intersect {1, ..., n}| efficiently
        last_n = max(memoized_intersection_size.keys())
        current_size = memoized_intersection_size[last_n]
        for i in range(last_n + 1, n + 1):
            if i in S_set:
                current_size += 1
            memoized_intersection_size[i] = current_size
        
        intersection_size = memoized_intersection_size[n]
        
        return cantor_pairing(n, intersection_size)

    return f_S

def main():
    """
    Main function to demonstrate the properties of f_S.
    """
    # Let S be the set of even numbers, as an example.
    # In a full proof, we would use a member of a large almost disjoint family.
    S = [i for i in range(2, 51, 2)]
    print(f"Using the sample set S = the first 25 even numbers: {S}")
    
    f_S = construct_f_S(S)
    
    print("\nVerifying that f_S(n) is strictly increasing for n = 1 to 20.")
    print("f_S(n) = pi(n, |S intersect {1, ..., n}|)\n")
    print("n | |S intersect {1..n}| |   pi(n, |...|) = f_S(n)")
    print("-" * 40)

    # Check for n from 1 up to a certain value
    limit = 20
    # Pre-calculate for n=0 to have a base case for comparison
    f_n_minus_1 = f_S(0)

    for n in range(1, limit + 1):
        # Calculate intersection size |S intersect {1, ..., n}|
        intersection_size = len(set(range(1, n + 1)).intersection(S))
        f_n = f_S(n)
        
        # We demonstrate the calculation: pi(n, intersection_size)
        k1 = n
        k2 = intersection_size
        result = (k1 + k2) * (k1 + k2 + 1) // 2 + k2
        
        print(f"{n:<2}| {k2:<21}| pi({k1}, {k2}) = {result}")

        # Assert that the function is strictly increasing
        if n > 0:
            assert f_n > f_n_minus_1, f"Error: f_S({n}) is not greater than f_S({n-1})"
        f_n_minus_1 = f_n

    print("\nVerification successful: The function f_S is strictly increasing.")
    print("This property is a key requirement for the functions used to construct the antichain.")
    print("The largest possible cardinality of such an antichain is 2^N_0.")


if __name__ == "__main__":
    main()
