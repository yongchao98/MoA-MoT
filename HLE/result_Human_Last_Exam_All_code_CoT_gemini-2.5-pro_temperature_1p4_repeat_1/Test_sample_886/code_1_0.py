import itertools

def solve_chess_problem():
    """
    Calculates the total number of distinct checkmate positions achievable
    by a single piece formed by combining the moves of two standard pieces.
    """

    # --- Step 1: Define base piece attack functions ---
    def can_attack_pawn(p_pos, t_pos):
        # Pawn attacks one step diagonally in any direction.
        return abs(p_pos[0] - t_pos[0]) == 1 and abs(p_pos[1] - t_pos[1]) == 1

    def can_attack_knight(p_pos, t_pos):
        dx = abs(p_pos[0] - t_pos[0])
        dy = abs(p_pos[1] - t_pos[1])
        return (dx == 1 and dy == 2) or (dx == 2 and dy == 1)

    def can_attack_bishop(p_pos, t_pos):
        if p_pos == t_pos: return False
        return abs(p_pos[0] - t_pos[0]) == abs(p_pos[1] - t_pos[1])

    def can_attack_rook(p_pos, t_pos):
        if p_pos == t_pos: return False
        return p_pos[0] == t_pos[0] or p_pos[1] == t_pos[1]

    def can_attack_queen(p_pos, t_pos):
        return can_attack_rook(p_pos, t_pos) or can_attack_bishop(p_pos, t_pos)

    def can_attack_king(p_pos, t_pos):
        return max(abs(p_pos[0] - t_pos[0]), abs(p_pos[1] - t_pos[1])) == 1

    base_pieces = {
        'P': can_attack_pawn, 'N': can_attack_knight, 'B': can_attack_bishop,
        'R': can_attack_rook, 'Q': can_attack_queen, 'K': can_attack_king,
    }

    # --- Step 2: Create the 15 piece combinations ---
    piece_names = ['P', 'N', 'B', 'R', 'Q', 'K']
    piece_pairs = list(itertools.combinations(piece_names, 2))

    total_checkmates = 0
    results = {}
    
    print("Calculating checkmate positions for each piece combination...")

    # --- Step 3: Main calculation loop ---
    for p1_name, p2_name in piece_pairs:
        p1_func = base_pieces[p1_name]
        p2_func = base_pieces[p2_name]
        
        # Create the combined attack function for the pair
        attack_func = lambda p, t: p1_func(p, t) or p2_func(p, t)
        
        pair_name = f"{p1_name}+{p2_name}"
        pair_mate_count = 0

        # King must be on a non-edge square for a 3x3 grid around it
        for kx in range(1, 7):
            for ky in range(1, 7):
                king_pos = (kx, ky)
                
                # Define the 3x3 block of squares to be attacked
                block_squares = []
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        block_squares.append((king_pos[0] + dx, king_pos[1] + dy))
                
                # Iterate through all possible positions for the attacking piece
                for px in range(8):
                    for py in range(8):
                        piece_pos = (px, py)
                        
                        # The piece cannot be within the 3x3 block
                        if piece_pos in block_squares:
                            continue
                            
                        # Check if the piece can attack all 9 squares
                        can_mate = True
                        for target_pos in block_squares:
                            if not attack_func(piece_pos, target_pos):
                                can_mate = False
                                break
                        
                        if can_mate:
                            pair_mate_count += 1
                            
        results[pair_name] = pair_mate_count
        total_checkmates += pair_mate_count

    # --- Step 4: Display the results ---
    print("\n--- Individual Results ---")
    equation_parts = []
    for pair_name, count in results.items():
        if count > 0:
            print(f"The number of checkmate positions for {pair_name} is: {count}")
            equation_parts.append(str(count))

    print("\n--- Final Calculation ---")
    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {total_checkmates}")
    
    print("\n--- Total ---")
    print(f"The total number of distinct checkmate positions is: {total_checkmates}")


if __name__ == '__main__':
    solve_chess_problem()