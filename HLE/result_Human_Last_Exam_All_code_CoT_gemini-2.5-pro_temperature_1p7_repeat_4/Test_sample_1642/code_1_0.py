import math

def solve_queen_problem():
    """
    Calculates the maximum number m of white and black queens that can coexist
    on a 16x16 board without attacking each other.
    """
    n = 16
    max_m = 0
    best_k = 0
    best_j = 0

    # Iterate through all possible row (k) and column (j) partitions
    for k in range(1, n):
        for j in range(1, n):
            # White queens are in the top-left k x j subgrid.
            # The square (0,0) is white. The number of white squares is ceil(area/2).
            white_area = k * j
            available_white = math.ceil(white_area / 2.0)

            # Black queens are in the bottom-right (16-k) x (16-j) subgrid.
            # This subgrid starts at (k, j).
            black_area = (n - k) * (n - j)
            
            # The color of the starting square (k,j) determines the counts.
            # If k+j is even, (k,j) is a white square. Number of black squares is floor(area/2).
            if (k + j) % 2 == 0:
                available_black = math.floor(black_area / 2.0)
            # If k+j is odd, (k,j) is a black square. Number of black squares is ceil(area/2).
            else:
                available_black = math.ceil(black_area / 2.0)

            current_m = min(available_white, available_black)
            if current_m > max_m:
                max_m = current_m
                best_k = k
                best_j = j
    
    print("The optimal strategy involves partitioning the board and separating queens by the color of their squares.")
    print(f"The best partition found is splitting the board by k={best_k} rows and j={best_j} columns.")
    print("\n--- Calculation for the Optimal Partition ---")

    # Recalculate for the best case to show the numbers
    white_area = best_k * best_j
    available_white = math.ceil(white_area / 2.0)
    
    black_area = (n - best_k) * (n - best_j)
    if (best_k + best_j) % 2 == 0:
        available_black = math.floor(black_area / 2.0)
    else:
        available_black = math.ceil(black_area / 2.0)
        
    print(f"White queens' subgrid is {best_k}x{best_j}, containing {int(available_white)} white squares.")
    print(f"Black queens' subgrid is {n-best_k}x{n-best_j}, containing {int(available_black)} black squares.")
    print(f"The maximum number m is the minimum of these two values.")
    print(f"m = min({int(available_white)}, {int(available_black)})")
    
    final_m = int(min(available_white, available_black))
    print(f"m = {final_m}")
    print(f"\nTherefore, the maximum number m is {final_m}.")

solve_queen_problem()