from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + + +",
    "+ X - @ - $ * + + +",
    "+ - @ X - - @ - + +",
    "+ - X - @ - - - + +",
    "+ - X - - @ - X + +",
    "+ - @ - - - - + + +",
    "+ - @ X - @ - - - +",
    "+ + X - $ - - - - +",
    "+ + + + + X - - - +",
    "+ + + + + + + + + +"
]

# Convert the initial state into a more manageable format
grid = [list(row.split()) for row in initial_state]
rows, cols = len(grid), len(grid[0])

# Find the initial position of the player
for r in range(rows):
    for c in range(cols):
        if grid[r][c] == '*':
            player_pos = (r, c)

# Define possible moves: (dr, dc) for (up, down, left, right)
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Function to check if a position is within the grid and not a wall
def is_valid(r, c):
    return 0 <= r < rows and 0 <= c < cols and grid[r][c] != '+'

# Function to check if a position is a goal
def is_goal(r, c):
    return grid[r][c] == 'X' or grid[r][c] == '$'

# Heuristic function: count the number of boxes not on goals
def heuristic(current_grid):
    return sum(1 for r in range(rows) for c in range(cols) if current_grid[r][c] == '@' and not is_goal(r, c))

# A* search to find the solution
def solve_sokoban():
    queue = []
    heappush(queue, (0, player_pos, grid, ""))
    visited = set()
    visited.add((player_pos, tuple(map(tuple, grid))))

    while queue:
        _, (pr, pc), current_grid, path = heappop(queue)

        # Check if all boxes are on goals
        if all(is_goal(r, c) for r in range(rows) for c in range(cols) if current_grid[r][c] == '@'):
            return path

        for move, (dr, dc) in moves.items():
            nr, nc = pr + dr, pc + dc

            if is_valid(nr, nc):
                # If the next position is a box, try to push it
                if current_grid[nr][nc] == '@' or current_grid[nr][nc] == '$':
                    nnr, nnc = nr + dr, nc + dc
                    if is_valid(nnr, nnc) and (current_grid[nnr][nnc] == '-' or current_grid[nnr][nnc] == 'X'):
                        # Create a new grid state
                        new_grid = [row[:] for row in current_grid]
                        # Move the player
                        new_grid[pr][pc] = '-' if current_grid[pr][pc] == '*' else 'X'
                        new_grid[nr][nc] = '*'
                        # Push the box
                        new_grid[nnr][nnc] = '@' if current_grid[nnr][nnc] == '-' else '$'
                        new_grid[nr][nc] = '-' if current_grid[nr][nc] == '@' else 'X'

                        new_state = ((nr, nc), tuple(map(tuple, new_grid)))
                        if new_state not in visited:
                            visited.add(new_state)
                            heappush(queue, (len(path) + 1 + heuristic(new_grid), (nr, nc), new_grid, path + move))
                # If the next position is empty or a goal, just move the player
                elif current_grid[nr][nc] == '-' or current_grid[nr][nc] == 'X':
                    new_grid = [row[:] for row in current_grid]
                    new_grid[pr][pc] = '-' if current_grid[pr][pc] == '*' else 'X'
                    new_grid[nr][nc] = '*'

                    new_state = ((nr, nc), tuple(map(tuple, new_grid)))
                    if new_state not in visited:
                        visited.add(new_state)
                        heappush(queue, (len(path) + 1 + heuristic(new_grid), (nr, nc), new_grid, path + move))

    return "No solution found"

# Find the solution
solution = solve_sokoban()
print(solution)