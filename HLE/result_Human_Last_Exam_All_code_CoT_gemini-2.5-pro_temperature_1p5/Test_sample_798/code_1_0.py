import math

def solve():
    """
    Finds the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572, based on Behrend's construction.
    """
    # The number of sets required only depends on the size of the interval,
    # specifically the largest number to be represented.
    # The start of the interval, 10001, is not relevant for this construction method.
    N_end = 42149572
    
    min_sets = float('inf')
    best_d = -1
    best_k = -1

    # Iterate through possible values of the base 'd'. A reasonable upper
    # bound is sufficient as the number of sets grows with d.
    for d in range(2, 200):
        # For a given d, calculate the minimum number of digits 'k' needed
        # for the Behrend set to be large enough to contain N_end.
        # The largest number in a Behrend set is (d-1) * ((2d)^k - 1) / (2d-1).
        # We need: (d-1) * ((2d)^k - 1) / (2d-1) >= N_end
        # This simplifies to: (2d)^k >= N_end * (2d-1)/(d-1) + 1
        if d == 1: continue # d must be >= 2
        
        # We need to solve for k
        val = N_end * (2 * d - 1) / (d - 1) + 1
        
        # k >= log(val) / log(2*d)
        k = math.ceil(math.log(val) / math.log(2 * d))

        # The number of sets in Behrend's partition is given by the
        # number of possible values for the sum of squares of the digits.
        # This sum can range from 0 to k*(d-1)^2.
        num_sets = k * (d - 1)**2 + 1

        if num_sets < min_sets:
            min_sets = num_sets
            best_d = d
            best_k = k

    print(f"The minimum number of sets is found using Behrend's construction.")
    print(f"The optimal parameters are d = {best_d} and k = {best_k}.")
    print("The final equation for the number of sets is:")
    print(f"Number of sets = k * (d-1)^2 + 1")
    print(f"Number of sets = {best_k} * ({best_d}-1)^2 + 1 = {min_sets}")
    print("\nResult:")
    print(min_sets)


solve()