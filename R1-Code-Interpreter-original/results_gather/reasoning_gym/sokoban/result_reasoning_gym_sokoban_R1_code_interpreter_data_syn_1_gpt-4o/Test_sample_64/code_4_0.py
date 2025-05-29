from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (2, 6),
    'boxes': {(2, 3), (2, 5), (3, 2), (4, 5), (5, 5), (6, 2)},
    'goals': {(1, 3), (2, 2), (2, 4), (3, 1), (6, 1), (6, 3)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
              (1, 0), (1, 1), (1, 6), (1, 7),
              (2, 0), (2, 1), (2, 6), (2, 7),
              (3, 0), (3, 6), (3, 7),
              (4, 0), (4, 6), (4, 7),
              (5, 0), (5, 6), (5, 7),
              (6, 0), (6, 6), (6, 7),
              (7, 0), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6), (7, 7)}
}

# Define the possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Function to check if a position is within bounds and not a wall
def is_valid_position(pos, walls):
    return pos not in walls

# Function to check if all boxes are on goals
def is_solved(boxes, goals):
    return boxes == goals

# BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state['player'], frozenset(initial_state['boxes']), '')])
    visited = set()

    while queue:
        player, boxes, path = queue.popleft()

        if is_solved(boxes, initial_state['goals']):
            return path

        if (player, boxes) in visited:
            continue

        visited.add((player, boxes))

        for move, (dx, dy) in moves.items():
            new_player = (player[0] + dx, player[1] + dy)

            if not is_valid_position(new_player, initial_state['walls']):
                continue

            if new_player in boxes:
                new_box = (new_player[0] + dx, new_player[1] + dy)

                if new_box in boxes or not is_valid_position(new_box, initial_state['walls']):
                    continue

                new_boxes = set(boxes)
                new_boxes.remove(new_player)
                new_boxes.add(new_box)
                queue.append((new_player, frozenset(new_boxes), path + move))
            else:
                queue.append((new_player, boxes, path + move))

    return None

solution = solve_sokoban(initial_state)
print(solution)