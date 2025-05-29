# Define the initial board state
board = [
    ['H', '.', 'B', 'B', 'C', 'C'],
    ['H', 'D', 'D', '.', '.', 'L'],
    ['I', 'A', 'A', '.', '.', 'L'],
    ['I', 'E', 'E', 'E', '.', 'L'],
    ['F', 'F', 'J', 'K', 'G', 'G'],
    ['.', '.', 'J', 'K', '.', '.']
]

# Define the moves
moves = [
    ('L', 1),  # Move L down
    ('E', -1), # Move E left
    ('D', -1), # Move D up
    ('B', -1), # Move B left
    ('C', -1), # Move C left
    ('A', 3)   # Move A right to the exit
]

# Function to apply moves to the board
def apply_moves(board, moves):
    for car, move in moves:
        # Find the car's position and orientation
        positions = [(r, c) for r in range(6) for c in range(6) if board[r][c] == car]
        if len(positions) > 1:
            # Determine orientation
            if positions[0][0] == positions[1][0]:  # Horizontal
                row = positions[0][0]
                if move > 0:  # Move right
                    for _ in range(move):
                        if positions[-1][1] + 1 < 6 and board[row][positions[-1][1] + 1] == '.':
                            for pos in reversed(positions):
                                board[pos[0]][pos[1] + 1] = car
                                board[pos[0]][pos[1]] = '.'
                else:  # Move left
                    for _ in range(-move):
                        if positions[0][1] - 1 >= 0 and board[row][positions[0][1] - 1] == '.':
                            for pos in positions:
                                board[pos[0]][pos[1] - 1] = car
                                board[pos[0]][pos[1]] = '.'
            else:  # Vertical
                col = positions[0][1]
                if move > 0:  # Move down
                    for _ in range(move):
                        if positions[-1][0] + 1 < 6 and board[positions[-1][0] + 1][col] == '.':
                            for pos in reversed(positions):
                                board[pos[0] + 1][pos[1]] = car
                                board[pos[0]][pos[1]] = '.'
                else:  # Move up
                    for _ in range(-move):
                        if positions[0][0] - 1 >= 0 and board[positions[0][0] - 1][col] == '.':
                            for pos in positions:
                                board[pos[0] - 1][pos[1]] = car
                                board[pos[0]][pos[1]] = '.'
    return board

# Apply the moves
final_board = apply_moves(board, moves)

# Print the final board state
for row in final_board:
    print(''.join(row))