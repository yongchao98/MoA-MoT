import math

def solve():
    """
    Calculates the minimum number of initially infected sites
    to infect an n x n grid based on the 3-neighbor rule.
    """
    n = 14
    
    # The reasoning is based on a perimeter argument.
    # Let I be the number of initially infected sites.
    # Let N be the total number of sites, N = n*n.
    # The initial perimeter P is at most 4*I.
    # The number of sites to infect is N - I.
    # Each infection reduces the perimeter by at least 2.
    # So, P >= 2 * (N - I).
    # Combining these gives 4*I >= 2*(N - I), which simplifies to 6*I >= 2*N or I >= N/3.
    
    total_sites = n * n
    
    # Derivation of the lower bound:
    # 4 * I >= 2 * (n^2 - I)
    # 4 * I >= 2 * n^2 - 2 * I
    # 6 * I >= 2 * n^2
    # I >= (2 * n^2) / 6
    # I >= n^2 / 3
    
    lower_bound_float = total_sites / 3
    
    # Since the number of sites must be an integer, we take the ceiling.
    min_sites = math.ceil(lower_bound_float)
    
    print(f"The grid size is n = {n}.")
    print(f"The total number of sites is n*n = {total_sites}.")
    print("The lower bound for the number of initial sites (I) is derived from the inequality 6*I >= 2*n*n.")
    print(f"This gives I >= n*n / 3.")
    print(f"I >= {total_sites} / 3")
    print(f"I >= {lower_bound_float}")
    print(f"As the number of sites must be an integer, the minimum number is ceil({lower_bound_float}) = {min_sites}.")
    print(f"Final Answer: The minimum number of initially infected sites is {min_sites}.")

solve()
