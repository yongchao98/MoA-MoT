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
            new_min_col = min_col + delta
            new_max_col = new_min_col + length - 1
            if new_min_col < 0 or new_max_col >= 6:
                return False
            # Check if the new position is empty
            if delta > 0:
                return self.board[row][new_max_col] == '.'
            else:
                return self.board[row][new_min_col] == '.'
        else:
            col = coords[0][1]
            min_row = coords[0][0]
            new_min_row = min_row + delta
            new_max_row = new_min_row + length - 1
            if new_min_row < 0 or new_max_row >= 6:
                return False
            # Check if the new position is empty
            if delta > 0:
                return self.board[new_max_row][col] == '.'
            else:
                return self.board[new_min_row][col] == '.'

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
            for _, col in coords:
                self.board[row][col + delta] = car
        else:
            col = coords[0][1]
            for row, _ in coords:
                self.board[row + delta][col] = car
        return True

    def is_solved(self):
        coords = self.get_car_coordinates('A')
        return coords and coords[-1][1] == 5  # Check if the rightmost part of A is at the exit

    def __str__(self):
        return '\n'.join([''.join(row) for row in self.board])

def solve_puzzle():
    initial_state = State()
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        board_str = str(current_state)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if current_state.is_solved():
            return moves
            
        # Try moving each car
        for car in 'ABCDEFGHI':  # Prioritize moving blocking cars
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