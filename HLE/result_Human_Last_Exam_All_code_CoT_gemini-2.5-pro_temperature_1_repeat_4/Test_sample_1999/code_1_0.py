def solve_spice_girls_chess_puzzle():
    """
    This script solves the Spice Girls chessboard puzzle by mapping their order
    from the 'Wannabe' rap bridge to the first rank of a chessboard.
    """

    # Step 1: Establish the order of the Spice Girls based on the rap bridge.
    # The order is interpreted as:
    # 1. "me" (the rapper) -> Mel B (Scary)
    # 2. "MC" -> Mel C (Sporty)
    # 3. "Em" -> Emma (Baby)
    # 4. "G" -> Geri (Ginger)
    # 5. "Easy V" -> Victoria (Posh)
    spice_girls_order = ["Scary", "Sporty", "Baby", "Ginger", "Posh"]

    # Step 2: Define the starting positions on the first rank of a chessboard.
    # The members start at the queenside Rook (a1) and continue along the rank.
    rank_files = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Step 3: Identify the target square.
    # The White Queen starts on the d-file.
    queen_start_file = 'd'

    print("The established order of the Spice Girls is:")
    print(spice_girls_order)
    print("\nMapping the members to the chessboard's first rank, starting from a1:")

    # The problem asks to show the relationship, which we can think of as an equation
    # where the position (index) maps to a member.
    placement = {}
    for i in range(len(spice_girls_order)):
        square = f"{rank_files[i]}1"
        member = f"{spice_girls_order[i]} Spice"
        placement[square] = spice_girls_order[i]
        # Outputting each "equation" or mapping step.
        print(f"Position {i + 1} ({square}) = {member}")

    # Step 4: Find the member on the Queen's starting square.
    queen_start_square = f"{queen_start_file}1"
    target_member_nickname = placement[queen_start_square]

    print(f"\nThe White Queen's starting square is {queen_start_square}.")
    print(f"The member standing on square {queen_start_square} is {target_member_nickname} Spice.")
    print("\nThe single word that precedes 'Spice' is:")
    print(target_member_nickname)

solve_spice_girls_chess_puzzle()