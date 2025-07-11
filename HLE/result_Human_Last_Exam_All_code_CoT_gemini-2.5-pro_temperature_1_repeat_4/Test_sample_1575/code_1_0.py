def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given adjacent and non-adjacent swap rules.
    """
    N = 100  # Total number of elements
    M = 5    # The step size for a free non-adjacent swap

    # The number of elements in each group
    num_per_group = N // M

    # Define the target group for each starting group (1-indexed)
    # An element starting in group `i` moves to a position `101-i`.
    # The target group is `(101-i) mod 5`.
    # Using 1-based indexing for groups:
    # i=1 -> (1-1)%5=0 -> group 5
    # i=2 -> (1-2)%5=-1 -> group 4
    # i=3 -> (1-3)%5=-2 -> group 3
    # i=4 -> (1-4)%5=-3 -> group 2
    # i=5 -> (1-5)%5=-4 -> group 1
    targets = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}

    # Helper function to calculate the shortest distance between two points on a circle
    def circle_dist(a, b, n):
        return min(abs(a - b), n - abs(a - b))

    # Calculate the sum of shortest distances for one element from each group
    dist_sum = 0
    dist_list = []
    for i in range(1, M + 1):
        dist = circle_dist(i, targets[i], M)
        dist_sum += dist
        dist_list.append(dist)

    # Calculate the total displacement for all N elements
    total_displacement = num_per_group * dist_sum

    # Each move (adjacent swap) accounts for a total displacement of 2
    total_moves = total_displacement // 2

    # --- Output the explanation and final equation ---
    print("The problem simplifies to moving 5 groups of elements between 5 positions on a circle.")
    print("The minimum number of moves is half the total displacement required.")
    print("\n--- Calculation Steps ---")
    
    # Build the strings for the equation
    dist_val_str = " + ".join(map(str, dist_list))
    
    print(f"1. Sum of unit distances for the 5 groups = {dist_val_str} = {dist_sum}")
    print(f"2. Total displacement = (Elements per group) * (Sum of unit distances)")
    print(f"   {num_per_group} * {dist_sum} = {total_displacement}")
    print(f"3. Minimum moves = (Total displacement) / 2")
    
    # Print the final equation with all the numbers
    print("\n--- Final Equation ---")
    print(f"Moves = ({num_per_group} * ({dist_val_str})) / 2")
    print(f"Moves = ({num_per_group} * {dist_sum}) / 2")
    print(f"Moves = {total_displacement} / 2")
    print(f"Moves = {total_moves}")


if __name__ == '__main__':
    solve_reversal_moves()