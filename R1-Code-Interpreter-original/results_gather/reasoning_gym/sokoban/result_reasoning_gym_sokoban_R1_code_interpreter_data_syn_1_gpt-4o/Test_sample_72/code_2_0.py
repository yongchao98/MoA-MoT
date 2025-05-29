from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = {
    'player': (3, 2),
    'boxes': {(2, 1), (2, 6), (3, 3), (4, 1), (4, 5), (5, 2), (6, 7)},
    'goals': {(1, 7), (2, 5), (4, 6), (5, 6), (6, 1)},
    'walls': {(0, i) for i in range(9)} | {(7, i) for i in range(9)} | {(i, 0) for i in range(8)} | {(i, 8) for i in range(8)}
}

# Define the possible moves
moves = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}

# Function to check if a move is valid
def is_valid_move(player, direction, boxes, walls):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    if new_player in walls:
        return False
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        if new_box in walls or new_box in boxes:
            return False
    return True

# Function to apply a move
def apply_move(state, move):
    player, boxes = state['player'], state['boxes']
    direction = moves[move]
    new_player = (player[0] + direction[0], player[1] + direction[1])
    new_boxes = set(boxes)
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    return {'player': new_player, 'boxes': new_boxes, 'goals': state['goals'], 'walls': state['walls']}

# Function to check if the puzzle is solved
def is_solved(state):
    return state['boxes'] == state['goals']

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(state):
    return sum(min(abs(box[0] - goal[0]) + abs(box[1] - goal[1]) for goal in state['goals']) for box in state['boxes'])

# A* search to find the solution
def solve_sokoban(initial_state):
    queue = [(heuristic(initial_state), 0, initial_state['player'], frozenset(initial_state['boxes']), "")]
    visited = set()
    while queue:
        _, cost, player, boxes, path = heappop(queue)
        if boxes == initial_state['goals']:
            return path
        if (player, boxes) in visited:
            continue
        visited.add((player, boxes))
        for move in moves:
            if is_valid_move(player, moves[move], boxes, initial_state['walls']):
                new_state = apply_move({'player': player, 'boxes': boxes, 'goals': initial_state['goals'], 'walls': initial_state['walls']}, move)
                new_cost = cost + 1
                heappush(queue, (new_cost + heuristic(new_state), new_cost, new_state['player'], frozenset(new_state['boxes']), path + move))
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)