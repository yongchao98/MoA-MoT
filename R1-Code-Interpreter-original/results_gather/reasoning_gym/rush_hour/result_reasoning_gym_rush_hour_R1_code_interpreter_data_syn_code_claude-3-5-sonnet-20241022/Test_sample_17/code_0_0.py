from collections import deque
import copy

# Define the board and vehicles
class State:
    def __init__(self):
        self.board = [
            ['B', 'B', 'B', 'F', '.', '.'],
            ['C', 'C', 'C', 'F', '.', '.'],
            ['.', 'A', 'A', 'F', '.', '.'],
            ['.', '.', 'E', '.', '.', '.'],
            ['.', '.', 'E', 'D', 'D', 'D'],
            ['.', '.', '.', '.', '.', '.']
        ]
        self.moves = []

    def __str__(self):
        return '\n'.join([''.join(row) for row in self.board])

    def get_hash(self):
        return str(self)

def can_move(board, row, col, is_vertical, length, delta):
    if is_vertical:
        new_row = row + delta
        if delta > 0:  # Moving down
            if new_row + length - 1 >= 6: return False
            return board[new_row + length - 1][col] == '.'
        else:  # Moving up
            if new_row < 0: return False
            return board[new_row][col] == '.'
    else:
        new_col = col + delta
        if delta > 0:  # Moving right
            if new_col + length - 1 >= 6: return False
            return board[row][new_col + length - 1] == '.'
        else:  # Moving left
            if new_col < 0: return False
            return board[row][new_col] == '.'

def move_vehicle(board, row, col, is_vertical, length, delta):
    vehicle = board[row][col]
    new_board = [row[:] for row in board]
    
    if is_vertical:
        for i in range(length):
            new_board[row + i][col] = '.'
        for i in range(length):
            new_board[row + delta + i][col] = vehicle
    else:
        for i in range(length):
            new_board[row][col + i] = '.'
        for i in range(length):
            new_board[row][col + delta + i] = vehicle
    
    return new_board

def find_vehicle(board, row, col):
    vehicle = board[row][col]
    if vehicle == '.': return None
    
    # Check horizontal
    length = 1
    while col + length < 6 and board[row][col + length] == vehicle:
        length += 1
    if length > 1:
        return (row, col, False, length)
    
    # Check vertical
    length = 1
    while row + length < 6 and board[row + length][col] == vehicle:
        length += 1
    if length > 1:
        return (row, col, True, length)
    
    return (row, col, False, 1)

def solve():
    initial = State()
    queue = deque([initial])
    seen = {initial.get_hash()}
    
    while queue:
        current = queue.popleft()
        
        # Check if solved
        if current.board[2][4] == 'A':
            return current.moves
        
        # Try all possible moves
        for row in range(6):
            for col in range(6):
                if current.board[row][col] != '.':
                    vehicle_info = find_vehicle(current.board, row, col)
                    if not vehicle_info: continue
                    
                    vrow, vcol, is_vertical, length = vehicle_info
                    if (row, col) != (vrow, vcol): continue
                    
                    # Try moving in both directions
                    for delta in [-1, 1]:
                        if can_move(current.board, row, col, is_vertical, length, delta):
                            new_state = State()
                            new_state.board = move_vehicle(current.board, row, col, is_vertical, length, delta)
                            new_state.moves = current.moves + [f"{current.board[row][col]}{'+' if delta > 0 else '-'}{abs(delta)}"]
                            
                            new_hash = new_state.get_hash()
                            if new_hash not in seen:
                                seen.add(new_hash)
                                queue.append(new_state)

solution = solve()
print(' '.join(solution))