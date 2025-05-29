def is_valid_move(start, end, board):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx, dy) in [(2, 1), (1, 2)] and board.get(end, ' ') == '.'

def swap_knights():
    board = {
        'A3': '.', 'B3': '.', 'C3': ' ', 'D3': '.',
        'A2': ' ', 'B2': 'w', 'C2': 'w', 'D2': '.',
        'A1': ' ', 'B1': '.', 'C1': 'B', 'D1': 'B'
    }
    
    target = {
        'B2': 'B', 'C2': 'B',
        'C1': 'w', 'D1': 'w'
    }
    
    moves = []
    knight_positions = {'w': ['B2', 'C2'], 'B': ['C1', 'D1']}
    move_order = ['w', 'B']
    
    def backtrack(turn):
        if all(board[pos] == target[pos] for pos in target):
            return True
        
        color = move_order[turn % 2]
        for i, start in enumerate(knight_positions[color]):
            for dx, dy in [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]:
                end = chr(ord(start[0]) + dx) + str(int(start[1]) + dy)
                if end in board and is_valid_move(start, end, board):
                    board[end] = color
                    board[start] = '.'
                    knight_positions[color][i] = end
                    moves.append(f"{color},{start},{end}")
                    
                    if backtrack(turn + 1):
                        return True
                    
                    # Backtrack
                    board[start] = color
                    board[end] = '.'
                    knight_positions[color][i] = start
                    moves.pop()
        
        return False
    
    if backtrack(0):
        return moves
    else:
        return "No"

print(swap_knights())