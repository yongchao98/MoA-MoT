from collections import deque

class Vehicle:
    def __init__(self, name, coords, is_horizontal):
        self.name = name
        self.coords = coords
        self.is_horizontal = is_horizontal
        self.length = len(coords)

class State:
    def __init__(self):
        # Initialize the 6x6 board
        self.board = [
            ['B', 'B', 'C', 'C', '.', 'x'],
            ['D', 'D', 'J', 'E', 'E', 'M'],
            ['.', '.', 'J', 'A', 'A', 'M'],
            ['x', '.', 'K', 'F', 'F', 'M'],
            ['.', '.', 'K', 'L', 'G', 'G'],
            ['.', 'H', 'H', 'L', 'I', 'I']
        ]
        self.vehicles = self.find_vehicles()

    def find_vehicles(self):
        vehicles = {}
        visited = set()
        
        for i in range(6):
            for j in range(6):
                if self.board[i][j].isalpha() and (i,j) not in visited:
                    name = self.board[i][j]
                    coords = [(i,j)]
                    visited.add((i,j))
                    
                    # Check horizontal
                    if j+1 < 6 and self.board[i][j+1] == name:
                        coords.append((i,j+1))
                        visited.add((i,j+1))
                        if j+2 < 6 and self.board[i][j+2] == name:
                            coords.append((i,j+2))
                            visited.add((i,j+2))
                        is_horizontal = True
                    # Check vertical
                    elif i+1 < 6 and self.board[i+1][j] == name:
                        coords.append((i+1,j))
                        visited.add((i+1,j))
                        if i+2 < 6 and self.board[i+2][j] == name:
                            coords.append((i+2,j))
                            visited.add((i+2,j))
                        is_horizontal = False
                    else:
                        is_horizontal = True
                    
                    vehicles[name] = Vehicle(name, coords, is_horizontal)
        
        return vehicles

    def can_move(self, vehicle, direction):
        if vehicle.is_horizontal:
            row = vehicle.coords[0][0]
            if direction > 0:  # right
                next_col = max(c[1] for c in vehicle.coords) + 1
                return next_col < 6 and self.board[row][next_col] == '.'
            else:  # left
                next_col = min(c[1] for c in vehicle.coords) - 1
                return next_col >= 0 and self.board[row][next_col] == '.'
        else:
            col = vehicle.coords[0][1]
            if direction > 0:  # down
                next_row = max(c[0] for c in vehicle.coords) + 1
                return next_row < 6 and self.board[next_row][col] == '.'
            else:  # up
                next_row = min(c[0] for c in vehicle.coords) - 1
                return next_row >= 0 and self.board[next_row][col] == '.'

    def move(self, vehicle_name, direction):
        new_state = State()
        new_state.board = [row[:] for row in self.board]
        vehicle = self.vehicles[vehicle_name]
        
        # Clear current position
        for i, j in vehicle.coords:
            new_state.board[i][j] = '.'
            
        # Set new position
        new_coords = []
        for i, j in vehicle.coords:
            if vehicle.is_horizontal:
                new_state.board[i][j + direction] = vehicle_name
                new_coords.append((i, j + direction))
            else:
                new_state.board[i + direction][j] = vehicle_name
                new_coords.append((i + direction, j))
                
        new_state.vehicles = new_state.find_vehicles()
        return new_state

    def is_solved(self):
        return self.board[2][4] == 'A'

    def __str__(self):
        return '\n'.join(''.join(row) for row in self.board)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

def solve():
    initial_state = State()
    visited = {str(initial_state)}
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state.is_solved():
            return moves
            
        for vehicle_name, vehicle in current_state.vehicles.items():
            for direction in [-1, 1]:
                if current_state.can_move(vehicle, direction):
                    next_state = current_state.move(vehicle_name, direction)
                    state_str = str(next_state)
                    
                    if state_str not in visited:
                        visited.add(state_str)
                        move = f"{vehicle_name}{'+' if direction > 0 else '-'}1"
                        queue.append((next_state, moves + [move]))
    
    return None

solution = solve()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")