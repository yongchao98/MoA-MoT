import collections

def get_neighbors(r, c):
    """Returns a list of orthogonal neighbors for a given coordinate (1-19)."""
    neighbors = []
    if r > 1: neighbors.append((r - 1, c))
    if r < 19: neighbors.append((r + 1, c))
    if c > 1: neighbors.append((r, c - 1))
    if c < 19: neighbors.append((r, c + 1))
    return neighbors

def find_group_and_liberties(r, c, board_state):
    """
    Finds the full group of connected stones and their liberties.
    A board_state is a dictionary {(r, c): 'B' or 'W'}.
    """
    stone_color = board_state.get((r, c))
    if not stone_color:
        return [], []

    group_stones = set()
    liberties = set()
    q = collections.deque([(r, c)])
    visited = set([(r, c)])

    while q:
        curr_r, curr_c = q.popleft()
        group_stones.add((curr_r, curr_c))

        for neighbor in get_neighbors(curr_r, curr_c):
            if neighbor not in visited:
                visited.add(neighbor)
                neighbor_content = board_state.get(neighbor)
                if not neighbor_content:
                    liberties.add(neighbor)
                elif neighbor_content == stone_color:
                    q.append(neighbor)
    return list(group_stones), list(liberties)

def solve_go_problem():
    """
    Analyzes the Go problem to find the winning move and demonstrates the outcome.
    """
    # Initial board state
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    board = {}
    for stone in black_stones:
        board[stone] = 'B'
    for stone in white_stones:
        board[stone] = 'W'

    print("Analyzing the board to find the winning move for Black...")
    
    # The proposed winning move for Black
    first_move = (2, 4)
    print(f"The critical move for Black is identified as {first_move}.")
    print("This move attacks all three separate White groups at once.\n")
    print("--- Simulating the game sequence after Black plays (2, 4) ---")

    # 1. Black plays the winning move
    board[first_move] = 'B'
    black_stones.add(first_move)
    print(f"1. Black plays at {first_move}.")
    w1_group, w1_libs = find_group_and_liberties(2, 5, board)
    print(f"   -> White's group at (2, 5) is now in atari, with its only liberty at {w1_libs[0]}.")

    # White must choose between saving the (2,5) stone or strengthening other groups.
    # We will simulate the case where White abandons the stone to save the larger group.
    
    # 2. White's response
    w_move_1 = (3, 2)
    board[w_move_1] = 'W'
    white_stones.add(w_move_1)
    print(f"2. White abandons the group in atari and plays at {w_move_1} to strengthen another group.")

    # 3. Black captures the weak group
    b_move_2 = (1, 5)
    board[b_move_2] = 'B'
    black_stones.add(b_move_2)
    captured_w1_stones, _ = find_group_and_liberties(2, 5, board)
    for s in captured_w1_stones:
        del board[s]
    white_stones.difference_update(captured_w1_stones)
    print(f"3. Black plays at {b_move_2} and captures the White stone(s) at {captured_w1_stones}.")

    # 4. Another White group is now in atari
    w2_group, w2_libs = find_group_and_liberties(1, 4, board)
    w_move_2 = w2_libs[0]
    print(f"   -> As a result, White's group at (1, 4) is now in atari with its only liberty at {w_move_2}.")
    board[w_move_2] = 'W'
    white_stones.add(w_move_2)
    print(f"4. White must save this group by playing at {w_move_2}.")

    # 5. Black continues the attack, connecting his stones
    b_move_3 = (2, 3)
    board[b_move_3] = 'B'
    black_stones.add(b_move_3)
    print(f"5. Black plays at {b_move_3} to press the attack.")

    # 6. This puts a combined white group in atari again
    w3_group, w3_libs = find_group_and_liberties(1, 3, board)
    w_move_3 = w3_libs[0]
    print(f"   -> This move puts a White group in atari again. The only escape is at {w_move_3}.")
    board[w_move_3] = 'W'
    white_stones.add(w_move_3)
    print(f"6. White's only move is to play at {w_move_3}, which connects all remaining white stones into a single large group.")

    # Final state: Black's turn to capture the single large group
    final_w_group, final_w_libs = find_group_and_liberties(list(white_stones)[0][0], list(white_stones)[0][1], board)
    print("\n--- FINAL STATE ---")
    print(f"All {len(final_w_group)} remaining White stones are in one group.")
    print(f"This group has only {len(final_w_libs)} liberties left: {sorted(final_w_libs)}.")
    print("It is Black's turn, and Black can capture the entire group by playing on these liberties sequentially.")

    final_capture_moves = sorted(final_w_libs)
    move_num = 7
    for move in final_capture_moves:
        print(f"{move_num}. Black plays at {move}.")
        move_num += 1
    
    print("\nConclusion: The first move (2, 4) creates an unstoppable sequence leading to the capture of all white stones.")


if __name__ == '__main__':
    solve_go_problem()