from collections import deque
import copy

# Define the board and cars
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

    def get_car_position(self, car):
        for i in range(6):
            for j in range(6):
                if self.board[i][j] == car:
                    return (i, j)
        return None

    def move_car(self, car, delta):
        pos = self.get_car_position(car)
        if not pos:
            return False
        
        new_board = copy.deepcopy(self.board)
        i, j = pos
        orientation = self.cars[car]['orientation']
        length = self.cars[car]['length']
        
        # Clear current position
        for k in range(length):
            if orientation == 'h':
                new_board[i][j+k] = '.'
            else:
                new_board[i+k][j] = '.'
        
        # Set new position
        if orientation == 'h':
            j += delta
            if j < 0 or j + length > 6:
                return False
            for k in range(length):
                if new_board[i][j+k] != '.':
                    return False
                new_board[i][j+k] = car
        else:
            i += delta
            if i < 0 or i + length > 6:
                return False
            for k in range(length):
                if new_board[i+k][j] != '.':
                    return False
                new_board[i+k][j] = car
                
        self.board = new_board
        return True

    def is_solved(self):
        return self.get_car_position('A')[1] == 4

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