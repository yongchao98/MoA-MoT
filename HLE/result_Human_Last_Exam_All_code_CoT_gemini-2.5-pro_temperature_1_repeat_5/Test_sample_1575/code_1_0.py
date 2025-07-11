import math

def solve_reversal_moves():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with special swap rules.
    """
    N = 100  # Total number of elements
    M = 5    # The modulo for position classes (from swap(i, i+5))

    if N % M != 0:
        print("Error: The number of elements must be a multiple of the class size.")
        return

    K = N // M  # Number of elements per class
    
    total_moves = 0
    move_components = []

    print("Calculating moves for each class of elements:")
    print("-" * 40)

    for start_class in range(M):
        # All elements starting in the same class go to the same destination class.
        # Let's find the destination for an element e_i where (i-1)%M = start_class.
        # Destination position is N+1-i. Destination class is ((N+1-i)-1)%M = (N-i)%M.
        # Since (i-1)%M = start_class, i = k*M + start_class + 1 for some k.
        # Dest class = (N - (k*M + start_class + 1)) % M
        # Dest class = (N - start_class - 1) % M
        # In Python, % can return negative, so we add M to ensure it's positive.
        dest_class = (N - start_class - 1 + M) % M
        
        # Calculate the shortest distance on the circle of classes
        dist = abs(start_class - dest_class)
        shortest_dist = min(dist, M - dist)
        
        moves_for_class = K * shortest_dist
        total_moves += moves_for_class
        
        # Store for final equation output
        move_components.append(f"{K} * {shortest_dist}")

        print(f"Elements starting in Class {start_class} (20 elements):")
        print(f"  - Destination Class: {dest_class}")
        print(f"  - Shortest distance: {shortest_dist} moves per element")
        print(f"  - Total moves for this group: {K} * {shortest_dist} = {moves_for_class}")
        print("-" * 40)

    equation_str = " + ".join(move_components)
    print("Final Calculation:")
    print(f"{equation_str} = {total_moves}")
    print(f"\nThe minimum number of moves required is: {total_moves}")


if __name__ == "__main__":
    solve_reversal_moves()