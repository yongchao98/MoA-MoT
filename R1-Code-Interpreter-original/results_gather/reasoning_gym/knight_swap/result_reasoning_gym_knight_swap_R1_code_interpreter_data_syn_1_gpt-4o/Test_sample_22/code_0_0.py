def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    positions = {
        'w1': 'A3', 'w2': 'D2',
        'B1': 'C2', 'B2': 'C1'
    }
    # Target positions
    targets = {
        'w1': 'C2', 'w2': 'C1',
        'B1': 'A3', 'B2': 'D2'
    }
    # Possible moves
    moves = []
    # Black moves first
    turn = 'B'
    # Simulate moves
    while positions != targets:
        moved = False
        for knight, pos in positions.items():
            if knight[0] == turn:
                for target in targets.values():
                    if target not in positions.values() and is_valid_move(pos, target):
                        moves.append(f"{knight},{pos},{target}")
                        positions[knight] = target
                        moved = True
                        break
            if moved:
                break
        # Alternate turn
        turn = 'w' if turn == 'B' else 'B'
        # If no move was made, it's impossible
        if not moved:
            return "No"
    return moves

# Find and print the moves
print(find_moves())