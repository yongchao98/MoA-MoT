import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle by analyzing the connectivity of the move graph.

    The method is to determine if the starting positions of black and white knights
    can be swapped. This is possible only if the distribution of each color of
    knight across the disconnected components of the board's move graph is identical.
    """

    # 1. Model the 4x3 board and knight moves. Squares are numbered 0-11.
    #    0  1  2
    #    3  4  5
    #    6  7  8
    #    9 10 11
    rows, cols = 4, 3
    moves = {
      0: [5, 7],    1: [6, 8],      2: [3, 7],
      3: [2, 10],   4: [9, 11],     5: [0, 10],
      6: [1, 8],    7: [0, 2, 6, 8], 8: [1, 7],
      9: [4],       10: [3, 5],     11: [4]
    }

    # 2. Find the connected components of the move graph using BFS.
    def find_components(nodes, edges):
        all_nodes = set(nodes)
        components = []
        while all_nodes:
            component = set()
            queue = collections.deque([all_nodes.pop()])
            component.add(queue[0])
            while queue:
                u = queue.popleft()
                for v in edges.get(u, []):
                    if v not in component:
                        component.add(v)
                        queue.append(v)
                        if v in all_nodes:
                            all_nodes.remove(v)
            components.append(frozenset(component))
        return components

    board_squares = list(range(rows * cols))
    connected_components = find_components(board_squares, moves)

    print("Analyzing the Knights Puzzle on a 4x3 board.")
    print(f"The board has {len(connected_components)} disconnected move regions (components):")
    for i, comp in enumerate(connected_components):
        print(f"  Component {i+1}: {sorted(list(comp))}")
    print("-" * 40)

    # 3. Define the five initial configurations.
    configs = {
        'A': {
            'black': {0, 3, 6, 9},
            'white': {2, 5, 8, 11}
        },
        'B': {
            'black': {1, 6, 8},
            'white': {4, 9, 11}
        },
        'C': {
            'black': {2, 5},
            'white': {0, 7}
        },
        'D': {
            'black': {4, 10},
            'white': {0, 6}
        },
        'E': {
            'black': {0, 3, 4},
            'white': {1, 2, 5}
        }
    }

    solvable_configs = []
    # 4. Analyze each configuration for solvability.
    for name, pos in configs.items():
        print(f"Analyzing Configuration {name}:")
        black_pos = pos['black']
        white_pos = pos['white']
        print(f"  Black knights start at: {sorted(list(black_pos))}")
        print(f"  White knights start at: {sorted(list(white_pos))}")
        
        # Calculate the distribution signature for each set of knights.
        black_sig = sorted([len(black_pos.intersection(c)) for c in connected_components])
        white_sig = sorted([len(white_pos.intersection(c)) for c in connected_components])
        
        print(f"  Black knights distribution over components: {black_sig}")
        print(f"  White knights distribution over components: {white_sig}")
        
        if black_sig == white_sig:
            print("  Result: The distributions match. The configuration is SOLVABLE.\n")
            solvable_configs.append(name)
        else:
            print("  Result: The distributions do not match. The configuration is UNSOLVABLE.\n")

    print("-" * 40)
    print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")


solve_knights_puzzle()