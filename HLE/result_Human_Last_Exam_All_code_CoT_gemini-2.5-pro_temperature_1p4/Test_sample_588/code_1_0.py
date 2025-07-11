def get_lexicographically_least_tuple():
    """
    Finds the lexicographically least tuple (a1, b1, ..., al, bl) that satisfies the problem conditions.
    """
    
    # Define a function to calculate the Euler characteristic
    def chi(a, b):
        return (2 - 2 * a) * (2 - 2 * b)

    # We deduced that the minimal length l is 2.
    # We need to find two pairs (a1, b1) and (a2, b2)
    # such that chi(a1, b1) + chi(a2, b2) = 0 and the resulting tuple is minimal.
    # The genera must be non-negative integers not equal to 1.
    
    limit = 5 # A search limit sufficient for finding the smallest values
    
    # Generate a list of candidate pairs (a,b) and their Euler characteristics.
    # A pair (a,b) is a candidate if a!=1 and b!=1.
    candidates = []
    for a in range(limit):
        if a == 1:
            continue
        for b in range(limit):
            if b == 1:
                continue
            # Only consider pairs where a <= b to avoid duplicates like (0,2) and (2,0) for now.
            if a <= b:
                pair_chi = chi(a, b)
                if pair_chi != 0:
                    candidates.append(((a,b), pair_chi))

    # Sort candidates based on the pair (a,b) to ensure we find the lexicographically smallest tuple
    candidates.sort(key=lambda x: x[0])
    
    # Find the first combination of two pairs whose chis sum to 0.
    best_tuple = None
    
    for i in range(len(candidates)):
        p1, chi1 = candidates[i]
        for j in range(i, len(candidates)):
            p2, chi2 = candidates[j]
            
            if chi1 + chi2 == 0:
                # We have a valid set of two manifolds.
                # The pairs are p1 and p2. We need to check both (p1, p2) and (p2, p1).
                # To get the lexicographically smallest tuple, we need to consider all combinations.
                # A more direct search is better.
                pass # This search strategy is a bit complex, let's simplify.
                
    # Direct search for the tuple (a1,b1,a2,b2)
    search_space = [i for i in range(limit) if i != 1]
    for a1 in search_space:
        for b1 in search_space:
            for a2 in search_space:
                for b2 in search_space:
                    # To ensure the final tuple is ordered lexicographically
                    t = tuple(sorted([(a1,b1), (a2,b2)]))
                    current_tuple = t[0] + t[1]
                    if current_tuple != (a1,b1,a2,b2):
                        continue

                    chi1 = chi(a1, b1)
                    chi2 = chi(a2, b2)
                    
                    if chi1 != 0 and chi1 + chi2 == 0:
                        # This is the first solution found, so it must be the lexicographically smallest.
                        a1_f, b1_f, a2_f, b2_f = current_tuple
                        
                        chi1_f = chi(a1_f, b1_f)
                        chi2_f = chi(a2_f, b2_f)
                        
                        print("The two manifolds are M(a1, b1) and M(a2, b2) where:")
                        print(f"(a1, b1) = ({a1_f}, {b1_f}) -> Euler characteristic = {chi1_f}")
                        print(f"(a2, b2) = ({a2_f}, {b2_f}) -> Euler characteristic = {chi2_f}")
                        
                        print("\nThe sum of their Euler characteristics is the condition for the connected sum to be full:")
                        print(f"(2-2*{a1_f})*(2-2*{b1_f}) + (2-2*{a2_f})*(2-2*{b2_f}) = {chi1_f} + ({chi2_f}) = {chi1_f + chi2_f}")

                        print("\nThe lexicographically least tuple is:")
                        # Print tuple as a string without spaces
                        print(str(current_tuple).replace(" ", ""))

                        return current_tuple

# Execute the function to get the answer
result = get_lexicographically_least_tuple()
