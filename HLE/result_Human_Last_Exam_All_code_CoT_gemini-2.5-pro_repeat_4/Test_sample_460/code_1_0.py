import math

def solve_suitable_subset_problem():
    """
    Calculates the smallest integer u based on a known theorem in combinatorial economics.
    
    The problem asks for the smallest integer u such that for any set of agents and preferences,
    a "suitable" subset of items O can be found. The parameters are:
    - m: number of items
    - t: a threshold for the number of agents assigned to an item
    - u: a threshold for the number of agents preferring an outside item

    The formula for the smallest u is given by u = C(m-1, k) * t,
    where k = ceil(m/2) - 1.
    """
    
    # Given parameters
    m = 4
    t = 20
    
    # The formula for the smallest u is derived from advanced results in fair division.
    # u = C(m-1, k) * t, where k = ceil(m/2) - 1.
    
    print("Calculating the smallest suitable u for m={} and t={}".format(m, t))
    print("-" * 50)
    
    # Step 1: Calculate k
    n_for_comb = m - 1
    k_for_comb = math.ceil(m / 2) - 1
    
    print("The formula is u = C(m-1, k) * t")
    print("First, we calculate the parameters for the combination function C(n, k):")
    print(f"n = m - 1 = {m} - 1 = {n_for_comb}")
    print(f"k = ceil(m/2) - 1 = ceil({m}/2) - 1 = {math.ceil(m/2)} - 1 = {k_for_comb}")
    
    # Step 2: Calculate the binomial coefficient C(n, k)
    comb_result = math.comb(n_for_comb, k_for_comb)
    print(f"\nNext, we calculate the binomial coefficient:")
    print(f"C({n_for_comb}, {k_for_comb}) = {comb_result}")
    
    # Step 3: Calculate the final value of u
    u = comb_result * t
    print(f"\nFinally, we calculate u:")
    print(f"u = C({n_for_comb}, {k_for_comb}) * t = {comb_result} * {t} = {u}")
    
    print("-" * 50)
    print(f"The smallest integer u is {u}.")

if __name__ == '__main__':
    solve_suitable_subset_problem()
