from collections import deque

# Initial positions
initial_state = (('w', 'C3'), ('w', 'D2'), ('B', 'B1'), ('B', 'D1'))
target_state = (('B', 'C3'), ('B', 'D2'), ('w', 'B1'), ('w', 'D1'))

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (ord(col) - ord('A'), int(row) - 1)

# Convert coordinates to position
def coord_to_pos(coord):
    col, row = coord
    return (chr(col + ord('A')), str(row + 1))

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# Generate possible moves for a knight
def generate_moves(knight, empty_squares):
    color, pos = knight
    x, y = pos_to_coord(pos)
    moves = []
    for dx, dy in knight_moves:
        new_coord = (x + dx, y + dy)
        if is_within_board(new_coord) and coord_to_pos(new_coord) in empty_squares:
            moves.append((color, pos, coord_to_pos(new_coord)))
    return moves

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()

        # Check if we reached the target state
        if set(current_state) == set(target_state):
            return moves

        # Find empty squares
        occupied = {pos for _, pos in current_state}
        empty_squares = {coord_to_pos((x, y)) for x in range(4) for y in range(3)} - occupied

        # Generate all possible moves
        for i, knight in enumerate(current_state):
            possible_moves = generate_moves(knight, empty_squares)
            for move in possible_moves:
                color, from_pos, to_pos = move
                new_state = list(current_state)
                new_state[i] = (color, to_pos)
                new_state = tuple(new_state)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, moves + [f"{color},{from_pos},{to_pos}"]))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)