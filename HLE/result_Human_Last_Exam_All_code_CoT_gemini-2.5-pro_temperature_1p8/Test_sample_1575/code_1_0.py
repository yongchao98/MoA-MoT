import math

def solve():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given allowed operations.
    """
    N = 100
    MOD = 5
    
    # The non-adjacent swap (i, i+5) partitions the 100 elements into 5 groups
    # based on their index modulo 5. Elements within a group can be rearranged for free.
    group_size = N // MOD
    
    # Reversing the sequence means an element at index i moves to 101 - i.
    # We analyze the required permutation of group contents.
    # An element from group j = i mod 5 must move to group k = (101 - i) mod 5.
    # This simplifies to the mapping j -> (1 - j) mod 5.
    # - Group 0 contents swap with Group 1 contents.
    # - Group 2 contents swap with Group 4 contents.
    # - Group 3 contents stay within Group 3 (cost = 0).
    
    # The groups are arranged cyclically. The cost to swap contents of two groups
    # depends on their distance 'd' in the cycle 0-1-2-3-4-0.
    # Cost = (2*d - 1) * group_size
    
    # 1. Cost to swap contents of Group 0 and Group 1
    # They are adjacent, so their distance d is 1.
    d_01 = 1
    cost_01 = (2 * d_01 - 1) * group_size
    
    # 2. Cost to swap contents of Group 2 and Group 4
    # The distance on the cycle is min(|2-4|, 5 - |2-4|) = 2.
    d_24 = 2
    cost_24 = (2 * d_24 - 1) * group_size

    # The total cost is the sum of costs for the independent swaps.
    total_moves = cost_01 + cost_24
    
    print(f"The sequence of {N} elements is partitioned into {MOD} groups of size {group_size}.")
    print("Reversing the sequence requires swapping the contents between certain groups.")
    print("The required swaps are (Group 0 <-> Group 1) and (Group 2 <-> Group 4).\n")
    
    print("Calculating cost for swapping contents of Group 0 and Group 1:")
    print(f"Distance on the cycle: d = {d_01}")
    print(f"Cost = (2 * {d_01} - 1) * {group_size} = {cost_01}\n")
    
    print("Calculating cost for swapping contents of Group 2 and Group 4:")
    print(f"Distance on the cycle: d = {d_24}")
    print(f"Cost = (2 * {d_24} - 1) * {group_size} = {cost_24}\n")

    print("The final equation for the total minimum moves is:")
    print(f"Total Moves = Cost(0,1) + Cost(2,4)")
    print(f"Total Moves = {cost_01} + {cost_24} = {total_moves}")

solve()