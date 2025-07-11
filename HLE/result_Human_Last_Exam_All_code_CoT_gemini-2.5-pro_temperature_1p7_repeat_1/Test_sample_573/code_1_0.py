import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given problem.
    """
    # Given parameters
    N = 200
    k = 13

    # The set of possible objects in the slice category is {k, k+1, ..., N}.
    # M is the size of this set.
    M = N - k + 1

    print(f"Given N = {N} and k = {k}.")
    print(f"The number of objects to choose from is M = {N} - {k} + 1 = {M}.")
    print("-" * 50)

    # List to store the results
    results = []

    # Calculate the number of n-simplices for n from 0 to 5.
    # The formula is the number of ways to choose a multiset of size n+1 from M items,
    # which is given by the binomial coefficient C(M + (n+1) - 1, n+1) = C(M+n, n+1).
    for n in range(6):
        # Number of simplices using the formula C(n, k) = n! / (k! * (n-k)!)
        # Here, the 'n' in C(n,k) is M+n, and 'k' is n+1.
        total_items = M + n
        items_to_choose = n + 1
        
        num_simplices = math.comb(total_items, items_to_choose)
        results.append(num_simplices)

        print(f"For n = {n}:")
        print(f"The number of {n}-simplices is given by C({M}+n, n+1) = C({total_items}, {items_to_choose}).")
        print(f"Result: {num_simplices}\n")
        
solve_simplices_count()