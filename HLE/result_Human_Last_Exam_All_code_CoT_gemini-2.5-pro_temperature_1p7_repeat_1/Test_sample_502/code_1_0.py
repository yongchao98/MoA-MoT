import math

def solve():
    """
    This function demonstrates the existence of a complex reflection group
    of order 10000 in GL(10, C) by constructing one as a product of
    cyclic groups, and thus finds the maximum possible dimension for R/I.
    """

    # The order of the group G
    group_order = 10000

    # The dimension of the vector space C^n
    n = 10

    # We can construct a reflection group as a direct product of n cyclic groups,
    # G = Z_k1 x Z_k2 x ... x Z_k10.
    # The order of this group is the product of the individual orders |G| = k1*k2*...*k10.
    # We need to find positive integers k_i whose product is 10000.
    # A simple choice is to have one k_i equal to 10000 and the others equal to 1.
    k_orders = [1] * n
    k_orders[0] = group_order

    # The resulting group order is the product of these values.
    # The python `math.prod` function calculates this product.
    calculated_order = math.prod(k_orders)

    # By the theorem mentioned, since we have found a complex reflection group
    # of order 10000, the maximum possible dimension of R/I is this order.
    max_dimension = calculated_order
    
    # We need to output the numbers in the final equation.
    # The equation is the product of the k_i values.
    equation_str = " * ".join(map(str, k_orders))
    
    print("The orders of the cyclic groups chosen are: {}".format(k_orders))
    print("The final equation for the group order is:")
    print("{} = {}".format(equation_str, max_dimension))
    print("\nThe largest possible dimension for R/I is {}".format(max_dimension))

solve()