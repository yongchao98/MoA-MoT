import math

def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    # The dimension of the vector space V corresponding to the defect group D of order 16=2^4.
    n = 4

    print("The problem is to find the maximum order of an odd-order subgroup E of GL(4, 2).")
    print("Let W be the fixed-point-free part of the vector space V under the action of E.")
    print("The order of E must divide 2^d - 1, where d is the dimension of W.")
    print("We check all possible dimensions for W from 1 to 4.")

    possible_orders = []
    
    # Iterate through possible dimensions of the subspace W
    for d in range(1, n + 1):
        # For each dimension d, the order of E must divide 2^d - 1
        divisor = (2**d) - 1
        possible_orders.append(divisor)
        print(f"If dim(W) = {d}, |E| must divide 2^{d} - 1 = {divisor}.")

    # The highest possible order is the maximum of these divisors.
    max_order = max(possible_orders)
    
    # Building the final equation string
    # E.g., max(1, 3, 7, 15) = 15
    equation_parts = ", ".join(map(str, sorted(possible_orders)))
    
    print("\nThe highest possible order for E is the maximum of these values.")
    print(f"Final calculation: max({equation_parts}) = {max_order}")
    print(f"\nThus, the highest order that E can have is {max_order}.")

solve()
