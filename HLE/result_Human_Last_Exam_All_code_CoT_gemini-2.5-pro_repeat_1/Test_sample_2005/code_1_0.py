import math

def solve_hyperknight_problem():
    """
    Calculates the minimum number of moves for the hyper-knight.
    """
    # Number of dimensions of the hypercube
    D = 7

    # Goal: Change each of the D coordinates from 0 to 2.
    # A change of +2 (mod 3) can be achieved by:
    #   - Two '+1' operations (2 total elementary ops)
    #   - One '-1' operation (1 total elementary op)

    # Let 'a' be the number of coordinates changed by -1 (1 op each).
    # Let 'b' be the number of coordinates changed by +2 (2 ops each).
    # We have a + b = D.

    # Total elementary operations = a * 1 + b * 2.
    # Each move consists of 2 elementary operations, so the total must be even.
    # a + 2b must be even. Since 2b is even, 'a' must be even.

    # Number of moves k = (a + 2b) / 2.
    # Substitute b = D - a: k = (a + 2(D - a)) / 2 = (2D - a) / 2.
    # To minimize k, we must maximize 'a'.

    # The maximum value for 'a' is the largest even integer <= D.
    if D % 2 == 0:
        max_a = D
    else:
        max_a = D - 1

    a = max_a
    b = D - a

    # Calculate the minimum number of moves
    min_moves = (a * 1 + b * 2) / 2

    # --- Output the reasoning and calculation ---
    print("The problem is to find the minimum number of moves to change 7 coordinates from 0 to 2.")
    print("A change of +2 (mod 3) can be achieved with one '-1' operation or two '+1' operations.")
    print("To minimize moves, we should maximize the use of the more efficient '-1' operation.")
    print("\nLet 'a' be the number of coordinates changed by '-1' (1 operation each).")
    print("Let 'b' be the number of coordinates changed by '+2' (2 operations each).")
    print(f"\nThe total number of operations is (1 * a + 2 * b). This must be even, so 'a' must be even.")
    print(f"To minimize moves, we must maximize 'a'. The maximum even number for 'a' where a <= {D} is {a}.")
    print(f"This gives us a = {a} and b = {D} - {a} = {b}.")
    
    print("\nThe final calculation for the minimum number of moves is:")
    print(f"Moves = (1 * a + 2 * b) / 2")
    print(f"Moves = (1 * {a} + 2 * {b}) / 2")
    print(f"Moves = ({a} + {2 * b}) / 2")
    print(f"Moves = {a + 2 * b} / 2")
    final_answer = int(min_moves)
    print(f"Moves = {final_answer}")
    return final_answer

if __name__ == '__main__':
    answer = solve_hyperknight_problem()
    # The final answer is wrapped according to the instruction.
    # print(f"\n<<<{answer}>>>")