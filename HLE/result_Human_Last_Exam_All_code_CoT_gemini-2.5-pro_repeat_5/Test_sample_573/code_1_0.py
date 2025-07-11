import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the specified category.
    """
    # Problem parameters
    N = 200
    k = 13

    # The objects of the overcategory are the integers i such that -N <= i <= k.
    # The number of these objects is M = k - (-N) + 1.
    M_base = N + k + 1

    print(f"Given N = {N} and k = {k}.")
    print("The category in question is the poset of integers from -N to k.")
    print(f"The number of objects in this category is M = {k} - (-{N}) + 1 = {M_base}.")
    print("\nAn n-simplex corresponds to a non-decreasing sequence of n+1 integers from this set.")
    print("The number of n-simplices is given by the multiset coefficient formula C(M + n, n + 1).")
    print("-" * 50)

    # Loop for n from 0 to 5
    for n in range(6):
        # The formula for the number of n-simplices is C(M_base + n, n + 1)
        upper_val = M_base + n
        lower_val = n + 1
        
        # Calculate the number of simplices using math.comb
        num_simplices = math.comb(upper_val, lower_val)
        
        # Print the result, showing the formula with numbers plugged in
        print(f"For n={n}, the number of {n}-simplices is:")
        print(f"C({N} + {k} + 1 + {n}, {n} + 1) = C({upper_val}, {lower_val}) = {num_simplices}")
        print()

if __name__ == "__main__":
    solve_simplices_count()