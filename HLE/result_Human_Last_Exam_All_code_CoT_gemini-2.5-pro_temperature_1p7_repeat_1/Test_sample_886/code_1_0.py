import collections

def solve_chess_puzzle():
    """
    Calculates the total number of distinct checkmate positions on a chessboard
    with a single piece that combines the moves of two different standard pieces.
    """
    Point = collections.namedtuple('Point', ['x', 'y'])
    ALL_SQUARES = [Point(x, y) for x in range(8) for y in range(8)]
    memo_moves = {}

    def get_moves(piece_char, pos):
        """
        A memoized function to get the set of squares a basic piece can attack from a given position.
        """
        if (piece_char, pos) in memo_moves:
            return memo_moves[(piece_char, pos)]

        moves = set()
        x, y = pos.x, pos.y

        if piece_char == 'K':
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < 8 and 0 <= ny < 8:
                        moves.add(Point(nx, ny))
        elif piece_char == 'R':
            for i in range(8):
                if i != y: moves.add(Point(x, i))
                if i != x: moves.add(Point(i, y))
        elif piece_char == 'B':
            for dx, dy in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
                for i in range(1, 8):
                    nx, ny = x + i * dx, y + i * dy
                    if 0 <= nx < 8 and 0 <= ny < 8:
                        moves.add(Point(nx, ny))
                    else:
                        break
        elif piece_char == 'N':
            deltas = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                      (2, 1), (2, -1), (-2, 1), (-2, -1)]
            for dx, dy in deltas:
                nx, ny = x + dx, y + dy
                if 0 <= nx < 8 and 0 <= ny < 8:
                    moves.add(Point(nx, ny))
        elif piece_char == 'P':
            # Pawn moves are directional. Assume 'forward' is towards increasing y-coordinate.
            # Combine one-step forward move and diagonal forward captures.
            if y + 1 < 8:
                moves.add(Point(x, y + 1))
                if x > 0: moves.add(Point(x - 1, y + 1))
                if x < 7: moves.add(Point(x + 1, y + 1))
        
        memo_moves[(piece_char, pos)] = moves
        return moves

    def get_combined_moves(piece_type, pos):
        """
        Gets the union of moves for a combined piece type (e.g., 'RN').
        """
        combined = set()
        for base_piece in piece_type:
            combined.update(get_moves(base_piece, pos))
        return combined

    mate_positions = set()
    base_pieces = ['R', 'B', 'N', 'K', 'P']
    
    # Generate the 10 combinations of two distinct piece types
    piece_combinations = []
    for i in range(len(base_pieces)):
        for j in range(i + 1, len(base_pieces)):
            piece_combinations.append(base_pieces[i] + base_pieces[j])

    for king_pos in ALL_SQUARES:
        king_escape_squares = get_moves('K', king_pos)
        for piece_pos in ALL_SQUARES:
            if king_pos == piece_pos:
                continue

            # King can capture an adjacent piece, which is not a checkmate.
            if piece_pos in king_escape_squares:
                continue

            for combo in piece_combinations:
                piece_attacks = get_combined_moves(combo, piece_pos)

                # Condition 1: The king must be in check.
                if king_pos not in piece_attacks:
                    continue

                # Condition 2: All king's escape squares must be attacked.
                if king_escape_squares.issubset(piece_attacks):
                    # This position is a mate. Add it to our set and check the next position.
                    position_tuple = (king_pos, piece_pos)
                    mate_positions.add(position_tuple)
                    break 

    print(f"{len(mate_positions)}")

solve_chess_puzzle()