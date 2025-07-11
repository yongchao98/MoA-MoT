def solve_prime_paths():
    """
    Calculates the number of distinct Prime Paths from (1, 1) to (5, 7) in 4 moves.

    The calculation is broken down by the number of horizontal (H) and vertical (V) moves.
    """
    
    print("Calculating the number of distinct Prime Paths...")
    print("-" * 30)

    # Case 1: 1 Horizontal Move, 3 Vertical Moves
    # This case has 4 arrangements (e.g., HVVV).
    # HVVV: 6 paths, VHVV: 9 paths, VVHV: 9 paths, VVVH: 6 paths.
    paths_1h_3v = 6 + 9 + 9 + 6
    print(f"Paths with 1 Horizontal and 3 Vertical moves: {paths_1h_3v}")

    # Case 2: 2 Horizontal Moves, 2 Vertical Moves
    # This case has C(4,2) = 6 arrangements (e.g., HVHV).
    # Each of the 6 arrangements (HHVV, HVHV, HVVH, VHHV, VHVH, VVHH) has 2*3=6 paths.
    paths_2h_2v = 6 * 6
    print(f"Paths with 2 Horizontal and 2 Vertical moves: {paths_2h_2v}")

    # Case 3: 3 Horizontal Moves, 1 Vertical Move
    # This case has 4 arrangements (e.g., HHHV).
    # Each arrangement has 2 paths.
    paths_3h_1v = 4 * 2
    print(f"Paths with 3 Horizontal and 1 Vertical move: {paths_3h_1v}")

    print("-" * 30)

    # Final Calculation
    total_paths = paths_1h_3v + paths_2h_2v + paths_3h_1v
    
    print("Total number of distinct Prime Paths is the sum of these cases.")
    print(f"Total = {paths_1h_3v} (from 1H,3V) + {paths_2h_2v} (from 2H,2V) + {paths_3h_1v} (from 3H,1V)")
    print(f"Total = {total_paths}")

solve_prime_paths()