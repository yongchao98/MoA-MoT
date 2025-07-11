import math

def solve_simplices():
    """
    Calculates the number of n-simplices for the given problem.
    """
    N = 200
    k = 13
    
    # The set of objects in the overcategory is {k, k+1, ..., N}
    # The size of this set is S.
    S = N - k + 1
    
    print(f"For N={N} and k={k}, the number of objects is S = {N} - {k} + 1 = {S}.")
    print("The number of n-simplices is given by the formula C(S+n, n+1).\n")

    for n in range(6):
        # We are choosing n+1 numbers from a set of size S with replacement.
        # This is equivalent to C(S + (n+1) - 1, n+1) = C(S+n, n+1).
        total_items = S + n
        items_to_choose = n + 1
        
        # Calculate the number of combinations
        result = math.comb(total_items, items_to_choose)
        
        # Build the equation string for clarity
        numerator_terms = [str(i) for i in range(total_items, total_items - items_to_choose, -1)]
        denominator_terms = [str(i) for i in range(items_to_choose, 0, -1)]
        
        num_str = " * ".join(numerator_terms)
        den_str = " * ".join(denominator_terms)
        
        if items_to_choose == 1:
             den_str = "1" # Cleaner display for C(n,1)
        
        print(f"For n = {n}, the number of {n}-simplices is C({total_items}, {items_to_choose}) = ({num_str}) / ({den_str}) = {result}")

solve_simplices()