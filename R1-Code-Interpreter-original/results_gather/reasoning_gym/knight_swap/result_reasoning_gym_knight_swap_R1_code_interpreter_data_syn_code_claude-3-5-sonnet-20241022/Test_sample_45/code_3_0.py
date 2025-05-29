from collections import deque
import copy

def print_board(board):
    for i, row in enumerate(board):
        print(f"{3-i} |", end=" ")
        print(" | ".join(cell if cell != '.' else ' ' for cell in row), end=" ")
        print("|")
        print("  " + "-" * 19)
    print("    A   B   C   D")

def make_move(board, move):
    color, from_pos, to_pos = move.split(',')
    from_col = ord(from_pos[0]) - ord('A')
    from_row = int(from_pos[1]) - 1
    to_col = ord(to_pos[0]) - ord('A')
    to_row = int(to_pos[1]) - 1
    
    new_board = copy.deepcopy(board)
    new_board[3-from_row][from_col] = '.'
    new_board[3-to_row][to_col] = color
    return new_board

def simulate_solution():
    board = [
        ['.', '.', 'B', '.'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    moves = ['w,B1,A3', 'B,C3,B1', 'w,D1,C3', 'B,B2,D1']
    
    print("Initial position:")
    print_board(board)
    
    for i, move in enumerate(moves, 1):
        board = make_move(board, move)
        print(f"\nAfter move {i}: {move}")
        print_board(board)
        
    # Check final positions
    white_positions = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'w')
    black_positions = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'B')
    
    print("\nFinal positions:")
    print(f"White knights at: {white_positions}")
    print(f"Black knights at: {black_positions}")
    
    return board

final_board = simulate_solution()