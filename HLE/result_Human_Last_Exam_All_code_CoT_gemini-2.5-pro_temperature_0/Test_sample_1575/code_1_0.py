import math

def solve():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    N = 100  # Total number of elements
    K = 5    # The distance of the free non-adjacent swap

    # Number of elements/positions in each group
    elements_per_group = N // K

    # The groups are indexed 0, 1, 2, 3, 4
    groups = list(range(K))
    total_moves = 0
    
    costs = []
    cost_expressions = []

    print("Calculating the cost for each group transition:")
    print("-" * 45)

    for g_initial in groups:
        # An element starting at position i (1-indexed) is in group (i-1) % K.
        # It must move to position (N+1 - i).
        # The destination group is ((N+1 - i) - 1) % K = (N - i) % K.
        # If (i-1)%K = g_initial, then i = m*K + g_initial + 1 for some integer m.
        # So, (N - i) % K = (N - (m*K + g_initial + 1)) % K
        #                  = (N - g_initial - 1) % K
        # For N=100, K=5, this is (99 - g_initial) % 5 = (4 - g_initial) % 5.
        g_final = (N - 1 - g_initial) % K

        # The distance is the minimum number of adjacent swaps to move between groups
        distance = min(abs(g_initial - g_final), K - abs(g_initial - g_final))

        # The cost for all elements in this group is elements_per_group * distance
        cost = elements_per_group * distance
        costs.append(cost)
        
        expression = f"{elements_per_group} * {distance}"
        cost_expressions.append(expression)

        print(f"Group {g_initial} -> Group {g_final}: Moves = {elements_per_group} * {distance} = {cost}")
        
        total_moves += cost

    print("-" * 45)
    
    # To meet the output requirement "output each number in the final equation"
    # we build the string carefully.
    sum_string = " + ".join(map(str, costs))
    print(f"Total moves = {sum_string} = {total_moves}")


solve()
<<<120>>>