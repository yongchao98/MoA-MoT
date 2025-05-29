from collections import deque

def sokoban_solver(initial_state):
    # Define the directions for movement
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # Convert the board into a mutable list of lists
    board = [list(row) for row in initial_state.split('\n')]
    
    # Find the initial position of the player and boxes
    player_pos = None
    boxes = set()
    goals = set()
    
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '*':
                player_pos = (i, j)
            elif cell == '@':
                boxes.add((i, j))
            elif cell == '$':
                boxes.add((i, j))
                goals.add((i, j))
            elif cell == 'X':
                goals.add((i, j))
    
    # Function to check if a position is within bounds and not a wall
    def is_valid(pos):
        x, y = pos
        return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'
    
    # Function to check if a move is valid
    def can_move(player, direction):
        dx, dy = directions[direction]
        new_player = (player[0] + dx, player[1] + dy)
        if not is_valid(new_player):
            return False
        if new_player not in boxes:
            return True
        new_box = (new_player[0] + dx, new_player[1] + dy)
        return is_valid(new_box) and new_box not in boxes
    
    # Function to apply a move
    def apply_move(player, direction):
        dx, dy = directions[direction]
        new_player = (player[0] + dx, player[1] + dy)
        if new_player in boxes:
            new_box = (new_player[0] + dx, new_player[1] + dy)
            boxes.remove(new_player)
            boxes.add(new_box)
        return new_player
    
    # Function to check if the puzzle is solved
    def is_solved():
        return all(box in goals for box in boxes)
    
    # BFS to find the solution
    queue = deque([(player_pos, '')])
    visited = set()
    visited.add((player_pos, frozenset(boxes)))
    
    while queue:
        player, path = queue.popleft()
        
        if is_solved():
            return path
        
        for direction in directions:
            if can_move(player, direction):
                new_player = apply_move(player, direction)
                new_state = (new_player, frozenset(boxes))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_player, path + direction))
                # Undo the move
                apply_move(new_player, 'UDLR'['UDLR'.index(direction) ^ 1])
    
    return "No solution found"

initial_state = """+ + + + + + + + + +  
+ * @ - - X X X - +  
+ @ $ @ - @ - @ - +  
+ X X - - - - @ - +  
+ + + @ - X + - - +  
+ + + X - - - - - +  
+ + + + + + + + + +  
+ + + + + + + + + +  """

solution = sokoban_solver(initial_state)
print(solution)