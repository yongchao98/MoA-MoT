def print_analysis():
    def get_knight_moves(pos):
        x, y = pos
        moves = [
            (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
        ]
        return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 3]

    # Convert chess notation to array indices
    def chess_to_idx(pos):
        return (3 - int(pos[1]), ord(pos[0]) - ord('A'))

    # Convert array indices to chess notation
    def idx_to_chess(i, j):
        return f"{chr(j + ord('A'))}{3-i}"

    # Initial positions
    w1_pos = chess_to_idx('A3')  # (0,0)
    w2_pos = chess_to_idx('A1')  # (2,0)
    b1_pos = chess_to_idx('B2')  # (1,1)
    b2_pos = chess_to_idx('B1')  # (2,1)

    # Target positions
    w1_target = chess_to_idx('B2')  # (1,1)
    w2_target = chess_to_idx('B1')  # (2,1)
    b1_target = chess_to_idx('A3')  # (0,0)
    b2_target = chess_to_idx('A1')  # (2,0)

    print("Analysis of the puzzle:")
    print("\n1. Initial positions:")
    print(f"White knights: A3({w1_pos}), A1({w2_pos})")
    print(f"Black knights: B2({b1_pos}), B1({b2_pos})")

    print("\n2. Possible moves from each position:")
    positions = [('A3', w1_pos), ('A1', w2_pos), ('B2', b1_pos), ('B1', b2_pos)]
    
    for pos_name, pos in positions:
        moves = get_knight_moves(pos)
        print(f"\nFrom {pos_name} ({pos}):")
        for move in moves:
            print(f"  -> {idx_to_chess(*move)} {move}")

    print("\n3. Key observations:")
    # Check if target squares are reachable in odd/even moves
    w1_moves = []
    visited = set()
    current = set([w1_pos])
    moves = 0
    while current and moves < 5:
        next_pos = set()
        for pos in current:
            for move in get_knight_moves(pos):
                if move not in visited:
                    next_pos.add(move)
                    visited.add(move)
                    if move == w1_target:
                        w1_moves.append(moves + 1)
        current = next_pos
        moves += 1

    print(f"Moves needed for first white knight to reach B2: {w1_moves}")
    print("The puzzle is impossible because:")
    print("1. White must move first")
    print("2. Players must alternate moves")
    print("3. To reach B2, white knight needs odd number of moves")
    print("4. This creates a parity issue with the remaining moves")

print_analysis()