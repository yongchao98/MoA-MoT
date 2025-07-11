import collections

# Constants for the board
BOARD_WIDTH = 3
BOARD_HEIGHT = 4
NUM_SQUARES = BOARD_WIDTH * BOARD_HEIGHT

# Global dictionary for pre-calculated moves
PRECALCULATED_MOVES = {}

def precompute_moves():
    """
    Pre-calculates all valid knight moves for each square on the board.
    The results are stored in the global PRECALCULATED_MOVES dictionary.
    """
    for r in range(BOARD_HEIGHT):
        for c in range(BOARD_WIDTH):
            pos = (r, c)
            PRECALCULATED_MOVES[pos] = []
            # All 8 possible knight move offsets
            offsets = [
                (1, 2), (1, -2), (-1, 2), (-1, -2),
                (2, 1), (2, -1), (-2, 1), (-2, -1)
            ]
            for dr, dc in offsets:
                new_r, new_c = r + dr, c + dc
                # Check if the new position is within the board boundaries
                if 0 <= new_r < BOARD_HEIGHT and 0 <= new_c < BOARD_WIDTH:
                    PRECALCULATED_MOVES[pos].append((new_r, new_c))

def is_solvable(initial_board):
    """
    Determines if a given Knights Puzzle configuration is solvable using Breadth-First Search (BFS).

    Args:
        initial_board (tuple): A tuple of length 12 representing the initial board state.
                               -1 for a black knight, 1 for a white knight, 0 for empty.

    Returns:
        bool: True if the configuration is solvable, False otherwise.
    """
    # The goal state is when all knights have swapped colors
    goal_board = tuple(-x for x in initial_board)

    # The state in our search is a tuple: (board_configuration, current_player_turn)
    # White moves first, represented by 1
    start_state = (initial_board, 1)
    
    # A queue for the BFS algorithm
    queue = collections.deque([start_state])
    # A set to keep track of visited states to avoid cycles and redundant work
    visited = {start_state}

    while queue:
        current_board, current_turn = queue.popleft()

        # Check if we have reached the goal configuration
        if current_board == goal_board:
            return True

        # Generate all possible next states from the current state
        for i, piece in enumerate(current_board):
            # Find a knight belonging to the player whose turn it is
            if piece == current_turn:
                r, c = i // BOARD_WIDTH, i % BOARD_WIDTH
                
                # Look up the pre-calculated valid moves for this knight's position
                for new_r, new_c in PRECALCULATED_MOVES.get((r, c), []):
                    new_i = new_r * BOARD_WIDTH + new_c
                    
                    # A knight can only move to an empty square
                    if current_board[new_i] == 0:
                        # Create the new board configuration after the move
                        new_board_list = list(current_board)
                        new_board_list[i], new_board_list[new_i] = 0, piece
                        new_board = tuple(new_board_list)
                        
                        # The turn alternates to the other player
                        next_turn = -current_turn
                        next_state = (new_board, next_turn)
                        
                        # If we haven't seen this state before, add it to the queue and visited set
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)
    
    # If the queue becomes empty and the goal was not found, the puzzle is unsolvable
    return False

def solve_puzzles():
    """
    Defines the five puzzle configurations from the image, solves them,
    and prints the result.
    """
    # Define board configurations based on the image
    # -1: Black Knight, 1: White Knight, 0: Empty Square
    # Board indices:
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    
    config_A = (-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1)

    board_B = [0] * NUM_SQUARES
    for i in [0, 6, 8]: board_B[i] = -1
    for i in [4, 9, 11]: board_B[i] = 1
    config_B = tuple(board_B)

    board_C = [0] * NUM_SQUARES
    for i in [2, 5]: board_C[i] = -1
    for i in [0, 7]: board_C[i] = 1
    config_C = tuple(board_C)

    board_D = [0] * NUM_SQUARES
    for i in [4, 10]: board_D[i] = -1
    for i in [0, 7]: board_D[i] = 1
    config_D = tuple(board_D)

    board_E = [0] * NUM_SQUARES
    for i in [0, 3, 4]: board_E[i] = -1
    for i in [1, 2, 5]: board_E[i] = 1
    config_E = tuple(board_E)

    configurations = {
        'A': config_A,
        'B': config_B,
        'C': config_C,
        'D': config_D,
        'E': config_E,
    }

    precompute_moves()

    solvable_configs = []
    for name, board in configurations.items():
        if is_solvable(board):
            solvable_configs.append(name)
    
    if solvable_configs:
        result = ", ".join(sorted(solvable_configs))
        print(f"The solvable configurations are: {result}")
    else:
        print("None of the configurations are solvable.")

# Execute the solver
solve_puzzles()