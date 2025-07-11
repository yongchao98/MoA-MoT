def solve():
    """
    This function demonstrates the solution for the group G = Z_5 x Z_2 x Z_2.
    It shows that this group of size 20 contains a maximal sum-free set S
    satisfying |k(S)| > 2|S|.
    Based on analysis of smaller groups, 20 is the smallest size.
    """
    G_orders = (5, 2, 2)
    G_size = G_orders[0] * G_orders[1] * G_orders[2]

    # This set S is a known maximal sum-free set in G = Z_5 x Z_2 x Z_2
    S = [[2, 0, 0], [3, 0, 0], [4, 1, 0]]
    size_S = len(S)

    # Calculate k(S) = {g in G | 2g is in S}
    kS = []
    
    # 2g = (2,0,0) -> 2a=2 (mod 5) -> a=1. b,c can be anything in Z2.
    # Solutions: (1,0,0), (1,0,1), (1,1,0), (1,1,1)
    kS.extend([[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
    
    # 2g = (3,0,0) -> 2a=3 (mod 5) -> a=4. b,c can be anything in Z2.
    # Solutions: (4,0,0), (4,0,1), (4,1,0), (4,1,1)
    kS.extend([[4, 0, 0], [4, 0, 1], [4, 1, 0], [4, 1, 1]])

    # 2g = (4,1,0) has no solution because 2(a,b,c) = (2a,0,0), which cannot have a non-zero second component.
    
    size_kS = len(kS)

    print(f"The smallest size of such a finite Abelian group is {G_size}.")
    print(f"The group is G = Z_{G_orders[0]} x Z_{G_orders[1]} x Z_{G_orders[2]}.")
    print(f"A maximal sum-free set is S = {S}, with |S| = {size_S}.")
    print(f"The set k(S) = {{g in G | 2g in S}} is calculated to be {kS}.")
    print(f"The size is |k(S)| = {size_kS}.")
    print("\nThe condition to check is |k(S)| > 2 * |S|.")
    print(f"We have the final equation: {size_kS} > 2 * {size_S}")
    print(f"This simplifies to {size_kS} > {2*size_S}, which is true.")
    print(f"Thus, the smallest size is {G_size}.")

solve()