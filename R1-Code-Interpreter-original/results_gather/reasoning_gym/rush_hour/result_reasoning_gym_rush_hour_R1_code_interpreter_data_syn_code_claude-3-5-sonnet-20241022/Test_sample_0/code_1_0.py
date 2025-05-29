from collections import deque
import copy

class State:
    def __init__(self):
        self.board = [
            ['.', '.', '.', 'G', '.', '.'],
            ['.', '.', '.', 'G', 'B', 'B'],
            ['.', '.', 'A', 'A', 'H', 'I'],
            ['x', '.', '.', '.', 'H', 'I'],
            ['F', 'C', 'C', 'D', 'D', 'I'],
            ['F', '.', '.', 'E', 'E', 'E']
        ]
        self.cars = {
            'A': {'orientation': 'h', 'length': 2},
            'B': {'orientation': 'h', 'length': 2},
            'C': {'orientation': 'h', 'length': 2},
            'D': {'orientation': 'v', 'length': 2},
            'E': {'orientation': 'h', 'length': 3},
            'F': {'orientation': 'v', 'length': 2},
            'G': {'orientation': 'v', 'length': 2},
            'H': {'orientation': 'v', 'length': 2},
            'I': {'orientation': 'v', 'length': 3}
        }

    def get_car_coordinates(self, car):
        coords = []
        for i in range(6):
            for j in range(6):
                if self.board[i][j] == car:
                    coords.append((i, j))
        return sorted(coords)

    def is_valid_move(self, car, delta):
        coords = self.get_car_coordinates(car)
        if not coords:
            return False
        
        orientation = self.cars[car]['orientation']
        length = self.cars[car]['length']
        
        if orientation == 'h':
            row = coords[0][0]
            min_col = coords[0][1]
            # Check new position bounds and obstacles
            new_min_col = min_col + delta
            if new_min_col < 0 or new_min_col + length > 6:
                return False
            # Check if new positions are empty
            for col in range(new_min_col, new_min_col + length):
                if col < min_col or col >= min_col + length:  # Only check new positions
                    if self.board[row][col] != '.':
                        return False
        else:  # vertical
            col = coords[0][1]
            min_row = coords[0][0]
            # Check new position bounds and obstacles
            new_min_row = min_row + delta
            if new_min_row < 0 or new_min_row + length > 6:
                return False
            # Check if new positions are empty
            for row in range(new_min_row, new_min_row + length):
                if row < min_row or row >= min_row + length:  # Only check new positions
                    if self.board[row][col] != '.':
                        return False
        return True

    def move_car(self, car, delta):
        if not self.is_valid_move(car, delta):
            return False
            
        coords = self.get_car_coordinates(car)
        orientation = self.cars[car]['orientation']
        
        # Clear current positions
        for i, j in coords:
            self.board[i][j] = '.'
            
        # Set new positions
        if orientation == 'h':
            row = coords[0][0]
            for idx, (_, col) in enumerate(coords):
                self.board[row][col + delta] = car
        else:
            col = coords[0][1]
            for idx, (row, _) in enumerate(coords):
                self.board[row + delta][col] = car
        return True

    def is_solved(self):
        coords = self.get_car_coordinates('A')
        return any(j == 4 for i, j in coords)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_state = State()
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        board_str = board_to_string(current_state.board)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if current_state.is_solved():
            return moves
            
        for car in current_state.cars:
            for delta in [-1, 1]:
                new_state = copy.deepcopy(current_state)
                if new_state.move_car(car, delta):
                    new_moves = moves + [f"{car}{'+' if delta > 0 else '-'}{abs(delta)}"]
                    queue.append((new_state, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")