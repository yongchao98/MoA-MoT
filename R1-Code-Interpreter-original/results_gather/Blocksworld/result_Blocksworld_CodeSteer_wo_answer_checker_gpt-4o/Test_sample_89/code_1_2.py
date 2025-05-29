def blocksworld_solution():
    # Initial state
    stacks = {
        1: ['G'],
        2: ['C', 'B'],
        3: ['E'],
        4: ['A', 'F', 'D']
    }

    # Goal state
    goal = {
        1: ['A', 'C'],
        2: ['B', 'D', 'E', 'G'],
        3: ['F'],
        4: []
    }

    # Moves to achieve the goal state
    moves = [
        ('D', 4, 2),
        ('F', 4, 3),
        ('A', 4, 1),
        ('E', 3, 2),
        ('G', 1, 2),
        ('C', 2, 1)
    ]

    # Print the moves
    for block, source, destination in moves:
        print(f"Move {block} from {source} to {destination}")

blocksworld_solution()