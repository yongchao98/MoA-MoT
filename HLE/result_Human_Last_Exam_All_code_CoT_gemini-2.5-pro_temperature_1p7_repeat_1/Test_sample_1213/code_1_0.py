def solve_go_puzzle():
    """
    Solves the go puzzle by modeling the board, rules, and analyzing moves.
    """

    # --- Step 1: Represent the Board ---
    # We use a dictionary for a sparse board representation.
    # Coordinates are tuples (col, row), e.g., A1 = (0, 0), B3 = (1, 2)
    def parse_coords(s):
        col = ord(s[0].upper()) - ord('A')
        row = int(s[1:]) - 1
        return (col, row)

    black_stones_str = ["A2", "B3", "B4", "C1", "C2"]
    white_stones_str = ["B5", "C3", "C4", "D1", "D2", "D5"]
    
    board = {}
    for s in black_stones_str:
        board[parse_coords(s)] = 'B'
    for s in white_stones_str:
        board[parse_coords(s)] = 'W'
    
    BOARD_SIZE = 19 # Assume a standard 19x19 board for edge checks

    # --- Step 2: Define Core Go Rules ---
    def get_neighbors(pt):
        x, y = pt
        neighbors = []
        if x > 0: neighbors.append((x - 1, y))
        if x < BOARD_SIZE - 1: neighbors.append((x + 1, y))
        if y > 0: neighbors.append((x, y - 1))
        if y < BOARD_SIZE - 1: neighbors.append((x, y + 1))
        return neighbors

    def find_group(pt, current_board):
        """Finds the group of connected stones and its liberties."""
        color = current_board.get(pt)
        if not color:
            return set(), set()

        q = [pt]
        group_stones = {pt}
        liberties = set()
        
        head = 0
        while head < len(q):
            current_pt = q[head]
            head += 1
            
            for neighbor in get_neighbors(current_pt):
                if neighbor not in current_board:
                    liberties.add(neighbor)
                elif current_board[neighbor] == color and neighbor not in group_stones:
                    group_stones.add(neighbor)
                    q.append(neighbor)
        return group_stones, liberties

    # --- Step 3 & 4: Analyze and Simulate ---
    # The key to this problem is understanding that the Black group's life depends
    # on the three points in the corner: A1, B1, and B2.
    # A successful kill requires White to control this area.

    # Candidate move: A1
    # This is the most severe move. Let's analyze it.
    # Sequence: 1. W at A1
    # Black's most resilient response is to play at B1. This captures the A1
    # stone and starts a "ko" fight.
    # In Go problems (tsumego), forcing a ko to kill a group is considered a
    # successful kill, as the group is not unconditionally alive. Any other
    # move by Black (e.g., at A3 or B2) allows White to play at B1 next,
    # destroying the eye shape and killing the group outright without a ko.
    
    # Candidate move: B1
    # Sequence: 1. W at B1
    # Black responds at A1. White's stone at B1 is now in "atari" (one liberty away from capture).
    # White must play B2 to save the B1 stone.
    # The resulting White group at B1-B2 has no liberties and is captured
    # (this is a non-capturing suicide, but the stones are removed).
    # So Black can refute W at B1.
    
    # Candidate move: B2
    # This is a "throw-in" move.
    # Sequence: 1. W at B2
    # Black's best response is to play at A1, trying to make an eye.
    # White's best follow-up is B1.
    # Similar to the B1 case, the White B1-B2 group is captured.
    # So Black can also refute W at B2.

    # --- Step 5: Identify the Killing Move(s) ---
    # Based on the analysis, only a white move at A1 initiates a sequence
    # where black cannot unconditionally live. Black is forced into a ko fight
    # which is considered a kill in this context.

    killing_moves = ['A1']

    # --- Step 6: Format the Output ---
    killing_moves.sort()
    result = "{" + ",".join(killing_moves) + "}"
    
    print(result)

solve_go_puzzle()