from collections import deque

class State:
    def __init__(self, board, moves=None):
        self.board = [list(row) for row in board]
        self.moves = moves if moves else []
        
    def get_car_positions(self):
        cars = {}
        for i in range(len(self.board)):
            for j in range(len(self.board[i])):
                if self.board[i][j] != '.' and self.board[i][j] != 'x':
                    if self.board[i][j] not in cars:
                        cars[self.board[i][j]] = [(i, j)]
                    else:
                        cars[self.board[i][j]].append((i, j))
        return cars

    def is_solved(self):
        cars = self.get_car_positions()
        if 'A' in cars:
            aa_pos = cars['A']
            return aa_pos[0][0] == 2 and max(pos[1] for pos in aa_pos) == len(self.board[0])-1
        return False

    def get_hash(self):
        return ''.join(''.join(row) for row in self.board)

    def is_horizontal(self, positions):
        return positions[0][0] == positions[1][0]

    def can_move(self, car_pos, direction):
        if self.is_horizontal(car_pos):
            row = car_pos[0][0]
            if direction > 0:
                col = max(pos[1] for pos in car_pos) + 1
                return col < len(self.board[0]) and self.board[row][col] == '.'
            else:
                col = min(pos[1] for pos in car_pos) - 1
                return col >= 0 and self.board[row][col] == '.'
        else:
            col = car_pos[0][1]
            if direction > 0:
                row = max(pos[0] for pos in car_pos) + 1
                return row < len(self.board) and self.board[row][col] == '.'
            else:
                row = min(pos[0] for pos in car_pos) - 1
                return row >= 0 and self.board[row][col] == '.'

    def move_car(self, car, car_pos, direction):
        new_board = [row[:] for row in self.board]
        for pos in car_pos:
            new_board[pos[0]][pos[1]] = '.'
        
        if self.is_horizontal(car_pos):
            row = car_pos[0][0]
            for pos in car_pos:
                new_board[row][pos[1] + direction] = car
        else:
            col = car_pos[0][1]
            for pos in car_pos:
                new_board[pos[0] + direction][col] = car
        
        return new_board

    def get_next_states(self):
        next_states = []
        cars = self.get_car_positions()
        
        # Prioritize moving cars that block the red car (AA)
        aa_pos = cars.get('A', [])
        if aa_pos:
            aa_row = aa_pos[0][0]
            blocking_cars = set()
            for col in range(max(pos[1] for pos in aa_pos) + 1, len(self.board[0])):
                if self.board[aa_row][col] not in ['.', 'x']:
                    blocking_cars.add(self.board[aa_row][col])
            
            # Try moving blocking cars first
            for car in blocking_cars:
                positions = cars[car]
                for direction in [-1, 1]:
                    if self.can_move(positions, direction):
                        new_board = self.move_car(car, positions, direction)
                        new_moves = self.moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        next_states.append(State(new_board, new_moves))

        # Then try moving other cars
        for car, positions in cars.items():
            if car not in blocking_cars:
                for direction in [-1, 1]:
                    if self.can_move(positions, direction):
                        new_board = self.move_car(car, positions, direction)
                        new_moves = self.moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        next_states.append(State(new_board, new_moves))
        
        return next_states

def solve_puzzle():
    initial_board = [
        ['x','.','G','B','B','J'],
        ['C','C','G','H','I','J'],
        ['F','A','A','H','I','K'],
        ['F','D','D','.','I','K'],
        ['E','E','.','x','.','.'],
        ['.','.','.','.','.','.']
    ]
    
    initial_state = State(initial_board)
    queue = deque([initial_state])
    visited = {initial_state.get_hash()}
    
    while queue:
        current_state = queue.popleft()
        
        if current_state.is_solved():
            return current_state.moves
            
        for next_state in current_state.get_next_states():
            state_hash = next_state.get_hash()
            if state_hash not in visited:
                visited.add(state_hash)
                queue.append(next_state)
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")