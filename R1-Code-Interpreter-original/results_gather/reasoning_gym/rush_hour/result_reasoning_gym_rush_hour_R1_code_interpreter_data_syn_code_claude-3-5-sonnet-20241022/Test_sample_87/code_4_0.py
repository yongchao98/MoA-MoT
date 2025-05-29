from collections import deque

class RushHourState:
    def __init__(self, board, moves=None):
        self.board = [list(row) for row in board]
        self.moves = moves if moves else []
        self.size = len(self.board)
        self.vehicles = self._get_vehicles()
        
    def _get_vehicles(self):
        vehicles = {}
        for i in range(self.size):
            for j in range(self.size):
                if self.board[i][j] != '.' and self.board[i][j] != 'x':
                    if self.board[i][j] not in vehicles:
                        vehicles[self.board[i][j]] = []
                    vehicles[self.board[i][j]].append((i, j))
        return {k: sorted(v) for k, v in vehicles.items()}
    
    def is_solved(self):
        red_car = self.vehicles['A']
        return max(x[1] for x in red_car) == self.size - 1
    
    def is_horizontal(self, vehicle):
        return self.vehicles[vehicle][0][0] == self.vehicles[vehicle][1][0]
    
    def can_move(self, vehicle, direction):
        positions = self.vehicles[vehicle]
        if self.is_horizontal(vehicle):
            row = positions[0][0]
            if direction > 0:  # right
                col = max(p[1] for p in positions) + 1
                return col < self.size and self.board[row][col] == '.'
            else:  # left
                col = min(p[1] for p in positions) - 1
                return col >= 0 and self.board[row][col] == '.'
        else:  # vertical
            col = positions[0][1]
            if direction > 0:  # down
                row = max(p[0] for p in positions) + 1
                return row < self.size and self.board[row][col] == '.'
            else:  # up
                row = min(p[0] for p in positions) - 1
                return row >= 0 and self.board[row][col] == '.'
    
    def move(self, vehicle, direction):
        new_board = [row[:] for row in self.board]
        positions = self.vehicles[vehicle]
        
        # Clear current positions
        for row, col in positions:
            new_board[row][col] = '.'
            
        # Set new positions
        is_horiz = self.is_horizontal(vehicle)
        for row, col in positions:
            new_row = row + (0 if is_horiz else direction)
            new_col = col + (direction if is_horiz else 0)
            new_board[new_row][new_col] = vehicle
            
        move_notation = f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"
        return RushHourState(new_board, self.moves + [move_notation])
    
    def get_blocking_score(self):
        red_car = self.vehicles['A']
        red_row = red_car[0][0]
        red_right = max(p[1] for p in red_car)
        score = 0
        for col in range(red_right + 1, self.size):
            if self.board[red_row][col] != '.':
                score += 1
        return score
    
    def get_hash(self):
        return ''.join(''.join(row) for row in self.board)

def solve_rush_hour():
    initial_board = [
        ['G', '.', '.', 'I', 'B', 'B'],
        ['G', '.', 'H', 'I', '.', 'J'],
        ['G', '.', 'H', 'A', 'A', 'J'],
        ['C', 'C', 'D', 'D', '.', 'K'],
        ['E', 'E', 'E', 'F', 'F', 'K'],
        ['.', '.', '.', '.', '.', '.']
    ]
    
    initial_state = RushHourState(initial_board)
    queue = deque([(initial_state, initial_state.get_blocking_score())])
    seen = {initial_state.get_hash()}
    
    while queue:
        current_state, _ = queue.popleft()
        
        if current_state.is_solved():
            return current_state.moves
        
        # Try all possible moves for each vehicle
        for vehicle in current_state.vehicles:
            for direction in [-1, 1]:
                if current_state.can_move(vehicle, direction):
                    new_state = current_state.move(vehicle, direction)
                    state_hash = new_state.get_hash()
                    
                    if state_hash not in seen:
                        seen.add(state_hash)
                        # Calculate priority based on blocking vehicles and move count
                        priority = new_state.get_blocking_score() + len(new_state.moves)
                        # Insert into queue based on priority
                        queue.append((new_state, priority))
        
        # Sort queue by priority
        queue = deque(sorted(queue, key=lambda x: x[1]))
    
    return None

solution = solve_rush_hour()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")