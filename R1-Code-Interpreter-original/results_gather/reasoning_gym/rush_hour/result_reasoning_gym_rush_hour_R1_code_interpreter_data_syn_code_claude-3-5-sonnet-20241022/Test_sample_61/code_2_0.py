from collections import deque, defaultdict

# Initial board state
initial_board = [
    list("BBBCCM"),
    list("DDEEEM"),
    list("AAJK.."),
    list("I.JKFF"),
    list("I.GGLx"),
    list("HH..L.")
]

class State:
    def __init__(self, board, parent=None, move=None):
        self.board = [row[:] for row in board]
        self.parent = parent
        self.move = move
        self.hash = ''.join(''.join(row) for row in board)
    
    def __eq__(self, other):
        return self.hash == other.hash
    
    def __hash__(self):
        return hash(self.hash)

def get_car_info(board):
    cars = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                cars[board[i][j]].append((i, j))
    
    car_orientations = {}
    for car, positions in cars.items():
        car_orientations[car] = 'H' if positions[0][0] == positions[-1][0] else 'V'
    
    return dict(cars), car_orientations

def get_blocking_cars(board, cars):
    if 'A' not in cars:
        return set()
    
    red_car = cars['A']
    red_row = red_car[0][0]
    red_right = max(pos[1] for pos in red_car)
    
    blocking = set()
    for col in range(red_right + 1, len(board[0])):
        if board[red_row][col] not in ['.', 'x']:
            blocking.add(board[red_row][col])
    return blocking

def get_moves(state, cars, orientations):
    moves = []
    blocking = get_blocking_cars(state.board, cars)
    
    for car, positions in cars.items():
        priority = 1 if car in blocking or car == 'A' else 2
        
        if orientations[car] == 'H':
            row = positions[0][0]
            left = min(pos[1] for pos in positions)
            right = max(pos[1] for pos in positions)
            
            # Try left
            if left > 0 and state.board[row][left-1] == '.':
                moves.append((priority, (car, -1)))
            # Try right
            if right < len(state.board[0])-1 and state.board[row][right+1] == '.':
                moves.append((priority, (car, 1)))
        else:
            col = positions[0][1]
            top = min(pos[0] for pos in positions)
            bottom = max(pos[0] for pos in positions)
            
            # Try up
            if top > 0 and state.board[top-1][col] == '.':
                moves.append((priority, (car, -1)))
            # Try down
            if bottom < len(state.board)-1 and state.board[bottom+1][col] == '.':
                moves.append((priority, (car, 1)))
    
    moves.sort()
    return [move[1] for move in moves]

def apply_move(state, cars, orientations, move):
    car, direction = move
    new_board = [row[:] for row in state.board]
    positions = cars[car]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    if orientations[car] == 'H':
        for pos in positions:
            new_board[pos[0]][pos[1] + direction] = car
    else:
        for pos in positions:
            new_board[pos[0] + direction][pos[1]] = car
    
    return State(new_board, state, move)

def is_solved(state, cars):
    if 'A' not in cars:
        return False
    red_car = cars['A']
    return max(pos[1] for pos in red_car) == 4

def solve():
    initial_state = State(initial_board)
    queue = deque([initial_state])
    visited = {initial_state}
    
    while queue:
        current = queue.popleft()
        cars, orientations = get_car_info(current.board)
        
        if is_solved(current, cars):
            moves = []
            while current.parent:
                moves.append(current.move)
                current = current.parent
            return moves[::-1]
        
        for move in get_moves(current, cars, orientations):
            next_state = apply_move(current, cars, orientations, move)
            if next_state not in visited:
                visited.add(next_state)
                queue.append(next_state)
    
    return None

# Solve and format the solution
solution = solve()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{car}{sign}1 "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")