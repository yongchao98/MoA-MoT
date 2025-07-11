import math

def get_minimum_perimeter(k):
    """
    Calculates the minimum possible perimeter for a polyomino of size k.
    The shape that minimizes perimeter is the one that is most "square-like".
    """
    if k == 0:
        return 0
    # The most compact shape is formed by arranging k squares in a grid
    # that is as close to a square as possible.
    # We can model this as a rectangle of width w and height h, with some
    # extra squares in a final row.
    w = int(math.sqrt(k))
    h = k // w
    remainder = k % w
    
    if remainder == 0:
        # Perfect rectangle of w x h
        perimeter = 2 * (w + h)
    else:
        # Rectangle of w x h with an extra partial row of 'remainder' squares.
        # This adds 2 to the perimeter of the w x h rectangle.
        perimeter = 2 * (w + h) + 2
        
    return perimeter

def solve_for_k():
    """
    Finds the maximal size K of the alien colony by checking the necessary
    condition P_min(K) <= 64 - K.
    """
    # We iterate downwards from the maximum possible squares to find the
    # largest K that satisfies the condition.
    for k in range(63, 0, -1):
        p_min = get_minimum_perimeter(k)
        v = 64 - k
        
        # The two main constraints for a valid final configuration are:
        # 1. p_min <= 30 (from max initial perimeter)
        # 2. p_min <= v (from the stopping condition)
        # We check both, but the second one is the binding constraint here.
        if p_min <= v and p_min <= 30:
            print(f"The maximal size of the colony (K) is {k}.")
            print("This is the largest K for which a shape can exist where the")
            print("minimum perimeter (P_min) is less than or equal to the number of vacant squares (V).")
            print("\nFinal equation check for K = " + str(k) + ":")
            print(f"P_min({k}) = {p_min}")
            print(f"V = 64 - {k} = {v}")
            # We must output each number in the final equation.
            # The key condition is P_min(K) <= V, or P_min(K) + K <= 64.
            # The maximal K will occur at the boundary where P_min(K) is close to V.
            final_sum = p_min + k
            print(f"The condition is P_min(K) + K <= 64.")
            print(f"Calculation: {p_min} + {k} = {final_sum}")
            return k

solve_for_k()