import collections

def solve_game():
    """
    Solves the game by generating the state graph and performing retrograde analysis.
    """
    
    # --- Constants and Piece Definitions ---
    BOARD_SIZE = 8
    P1_PIECES = frozenset(['K1', 'N1', 'R1'])
    P2_PIECES = frozenset(['K2', 'N2', 'R2'])

    # --- State Representation Helpers ---
    def board_dict_to_tuple(board_dict):
        return tuple(sorted(board_dict.items()))

    def board_tuple_to_dict(board_tuple):
        return dict(board_tuple)

    def board_dict_to_list(board_dict):
        board = [' '] * BOARD_SIZE
        for piece, pos in board_dict.items():
            if pos is not None:
                board[pos] = piece
        return board

    # --- Game Logic ---
    def get_king_moves(pos):
        moves = []
        for d in [-1, 1]:
            if 0 <= pos + d < BOARD_SIZE:
                moves.append(pos + d)
        return moves

    def get_knight_moves(pos):
        moves = []
        for d in [-2, 2]:
            if 0 <= pos + d < BOARD_SIZE:
                moves.append(pos + d)
        return moves

    def get_rook_moves(board_list, pos):
        moves = []
        # Right
        for new_pos in range(pos + 1, BOARD_SIZE):
            moves.append(new_pos)
            if board_list[new_pos] != ' ': break
        # Left
        for new_pos in range(pos - 1, -1, -1):
            moves.append(new_pos)
            if board_list[new_pos] != ' ': break
        return moves

    def is_in_check(board_list, player):
        king_piece = 'K1' if player == 1 else 'K2'
        rook_piece = 'R2' if player == 1 else 'R1'

        try:
            king_pos = board_list.index(king_piece)
            rook_pos = board_list.index(rook_piece)
        except ValueError:
            return False

        start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
        path_clear = all(board_list[i] == ' ' for i in range(start + 1, end))
        return path_clear

    def get_legal_next_states(board_tuple, player):
        board_dict = board_tuple_to_dict(board_tuple)
        board_list = board_dict_to_list(board_dict)
        player_pieces_set = P1_PIECES if player == 1 else P2_PIECES
        
        next_states = []
        
        pieces_to_move = [item for item in board_dict.items() if item[0] in player_pieces_set]

        for piece, pos in pieces_to_move:
            piece_type = piece[0]
            
            possible_dests = []
            if piece_type == 'K': possible_dests = get_king_moves(pos)
            elif piece_type == 'N': possible_dests = get_knight_moves(pos)
            elif piece_type == 'R': possible_dests = get_rook_moves(board_list, pos)
            
            for dest in possible_dests:
                dest_piece = board_list[dest]
                if dest_piece in player_pieces_set:
                    continue
                
                next_board_dict = board_dict.copy()
                next_board_dict[piece] = dest
                if dest_piece != ' ':
                    del next_board_dict[dest_piece]
                    
                next_board_list = board_dict_to_list(next_board_dict)

                if not is_in_check(next_board_list, player):
                    next_board_tuple = board_dict_to_tuple(next_board_dict)
                    next_player = 3 - player
                    next_states.append((next_board_tuple, next_player))
        return next_states

    # --- Main Solver Logic ---
    # 1. Generate state graph (all reachable states)
    initial_board_dict = {'K1': 0, 'N1': 1, 'R1': 2, 'R2': 5, 'N2': 6, 'K2': 7}
    initial_state = (board_dict_to_tuple(initial_board_dict), 1)

    q = collections.deque([initial_state])
    visited = {initial_state}
    predecessors = collections.defaultdict(list)
    all_states = [initial_state]
    move_counts = {}

    head = 0
    while head < len(all_states):
        state = all_states[head]
        head += 1
        board_tuple, player = state

        next_states = get_legal_next_states(board_tuple, player)
        move_counts[state] = len(next_states)
        
        for next_state in next_states:
            predecessors[next_state].append(state)
            if next_state not in visited:
                visited.add(next_state)
                all_states.append(next_state)

    # 2. Retrograde Analysis
    outcomes = {}  # state -> (outcome, depth in half-moves)
    analysis_q = collections.deque()
    loss_depth_precursors = collections.defaultdict(list)

    for state in all_states:
        if move_counts[state] == 0:
            outcomes[state] = ('LOSS', 0)
            analysis_q.append(state)

    while analysis_q:
        state = analysis_q.popleft()
        outcome, depth = outcomes[state]

        for pred_state in predecessors[state]:
            if pred_state in outcomes:
                continue
            
            if outcome == 'LOSS':
                outcomes[pred_state] = ('WIN', depth + 1)
                analysis_q.append(pred_state)
            else: # outcome == 'WIN'
                loss_depth_precursors[pred_state].append(depth)
                if len(loss_depth_precursors[pred_state]) == move_counts[pred_state]:
                    max_opp_win_depth = max(loss_depth_precursors[pred_state])
                    outcomes[pred_state] = ('LOSS', max_opp_win_depth + 1)
                    analysis_q.append(pred_state)

    # 3. Get and print the result
    final_outcome, final_depth = outcomes.get(initial_state, ('NO_FORCED_WIN', -1))

    if final_outcome == 'WIN':
        turns = (final_depth + 1) // 2
        print(turns)
    else:
        print("Player 1 cannot force a win against optimal play.")

solve_game()