import collections

def solve_linear_chess():
    """
    Solves a custom 1D chess-like game to find the shortest forced win for Player 1.
    The method uses retrograde analysis on the game's state space.
    """

    # --- Game Rules and Constants ---
    INITIAL_BOARD = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    PIECES = {'K1': 1, 'N1': 1, 'R1': 1, 'K2': 2, 'N2': 2, 'R2': 2}
    PLAYERS = {1: ('K1', 'N1', 'R1'), 2: ('K2', 'N2', 'R2')}
    BOARD_SIZE = 8

    # --- Memoization Caches ---
    memo_is_check = {}
    memo_legal_moves = {}

    # --- Helper Functions ---
    def is_check(board, player):
        """Checks if the given player's king is in check by an opponent's rook."""
        state = (board, player)
        if state in memo_is_check:
            return memo_is_check[state]

        player_king = f'K{player}'
        opponent_rook = f'R{3-player}'

        try:
            king_pos = board.index(player_king)
        except ValueError:
            memo_is_check[state] = False
            return False

        try:
            rook_pos = board.index(opponent_rook)
        except ValueError:
            memo_is_check[state] = False
            return False

        start, end = sorted((king_pos, rook_pos))
        for i in range(start + 1, end):
            if board[i] != ' ':
                memo_is_check[state] = False
                return False
        
        memo_is_check[state] = True
        return True

    def get_legal_moves(board, player):
        """Generates all legal moves for a given player."""
        state = (board, player)
        if state in memo_legal_moves:
            return memo_legal_moves[state]

        legal_moves = []
        player_pieces = PLAYERS[player]
        
        for pos, piece in enumerate(board):
            if piece not in player_pieces:
                continue

            if piece.startswith('K'):
                moves = [-1, 1]
            elif piece.startswith('N'):
                moves = [-2, 2]
            else: # Rook
                moves = []
            
            if piece.startswith('K') or piece.startswith('N'):
                for move in moves:
                    new_pos = pos + move
                    if 0 <= new_pos < BOARD_SIZE and (board[new_pos] == ' ' or PIECES.get(board[new_pos]) != player):
                        new_board_list = list(board)
                        new_board_list[pos], new_board_list[new_pos] = ' ', piece
                        new_board = tuple(new_board_list)
                        if not is_check(new_board, player):
                            legal_moves.append(new_board)
            elif piece.startswith('R'):
                for direction in [-1, 1]:
                    for i in range(1, BOARD_SIZE):
                        new_pos = pos + i * direction
                        if not (0 <= new_pos < BOARD_SIZE):
                            break
                        
                        target_cell = board[new_pos]
                        if target_cell == ' ':
                            new_board_list = list(board)
                            new_board_list[pos], new_board_list[new_pos] = ' ', piece
                            new_board = tuple(new_board_list)
                            if not is_check(new_board, player):
                                legal_moves.append(new_board)
                        elif PIECES.get(target_cell) != player: # Capture
                            new_board_list = list(board)
                            new_board_list[pos], new_board_list[new_pos] = ' ', piece
                            new_board = tuple(new_board_list)
                            if not is_check(new_board, player):
                                legal_moves.append(new_board)
                            break
                        else: # Friendly piece
                            break
        
        unique_moves = list(set(legal_moves))
        memo_legal_moves[state] = unique_moves
        return unique_moves

    # 1. Explore the entire reachable state space
    moves_from = collections.defaultdict(list)
    p1_turn_boards = set()
    p2_turn_boards = set()
    queue_for_graph = collections.deque([(INITIAL_BOARD, 1)])
    all_states = {(INITIAL_BOARD, 1)}

    while queue_for_graph:
        board, player = queue_for_graph.popleft()
        
        if player == 1:
            p1_turn_boards.add(board)
        else:
            p2_turn_boards.add(board)
            
        next_player = 3 - player
        legal_next_boards = get_legal_moves(board, player)
        moves_from[board] = legal_next_boards
        
        for next_board in legal_next_boards:
            if (next_board, next_player) not in all_states:
                all_states.add((next_board, next_player))
                queue_for_graph.append((next_board, next_player))

    # 2. Retrograde Analysis
    W1 = set() # Boards where P1 (to move) has a forced win.
    L2 = set() # Boards where P2 (to move) has a forced loss.
    
    for k in range(1, len(all_states) + 1):
        W1_new = set()
        L2_new = set()

        # Find P1 positions that are now winning (can move to a known P2 loss or checkmate).
        for board in p1_turn_boards:
            if board in W1:
                continue
            for next_board in moves_from[board]:
                if 'K2' not in next_board or next_board in L2:
                    W1_new.add(board)
                    break
        
        if INITIAL_BOARD in W1_new:
            game_length = 2 * k - 1
            print(game_length)
            return

        # Find P2 positions that are now losing (all moves go to a known P1 win).
        all_p1_wins = W1.union(W1_new)
        for board in p2_turn_boards:
            if board in L2:
                continue
            
            next_boards = moves_from[board]
            if next_boards and all(nb in all_p1_wins for nb in next_boards):
                L2_new.add(board)

        if not W1_new and not L2_new:
            break
            
        W1.update(W1_new)
        L2.update(L2_new)

    print("No forced win found for Player 1.")

solve_linear_chess()