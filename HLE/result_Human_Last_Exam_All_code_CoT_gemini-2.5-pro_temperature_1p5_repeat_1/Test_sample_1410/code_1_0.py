import chess
import chess.engine
import sys
import asyncio

async def find_mate_in_n():
    """
    Analyzes a chess FEN string to find the shortest mate sequence.
    """
    # The user might need to install stockfish, e.g., 'apt-get install stockfish' on Debian/Ubuntu
    # or 'brew install stockfish' on macOS.
    # On Windows, the executable needs to be downloaded and its path provided.
    engine_path = "stockfish"
    if sys.platform == "win32":
        engine_path += ".exe"

    try:
        # For Windows, you might need the full path, e.g., r"C:\path\to\stockfish.exe"
        transport, engine = await chess.engine.popen_uci(engine_path)
    except FileNotFoundError:
        print(f"Error: The '{engine_path}' engine could not be found.")
        print("Please install Stockfish and ensure it's in your system's PATH,")
        print("or specify the full path to the executable in the script.")
        # Gracefully exit if the engine is not available.
        await engine.quit()
        return

    # Set up the board from the FEN provided by the user
    board = chess.Board("4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20")

    # Analyze the position for a short time to find the forced mate
    try:
        result = await engine.play(board, chess.engine.Limit(time=2.0))
    except chess.engine.EngineTerminatedError:
        print("Stockfish engine terminated unexpectedly. Cannot complete analysis.")
        await engine.quit()
        return

    # A mate is guaranteed to be found if one exists
    if result.move:
        # Find the full line until mate
        pv = [result.move]
        board.push(result.move)
        while not board.is_game_over():
            try:
                analysis = await engine.analyse(board, chess.engine.Limit(depth=20))
                # Check if a mate sequence is found
                if 'pv' in analysis and analysis['pv']:
                    next_move = analysis['pv'][0]
                    pv.append(next_move)
                    board.push(next_move)
                else:
                    break # Stop if no clear best move is found
            except chess.engine.EngineTerminatedError:
                print("Stockfish engine terminated unexpectedly during analysis.")
                break


    # Reset board to print the sequence with standard algebraic notation
    board.reset()
    
    # Print the full winning line
    move_list_san = []
    move_num = board.fullmove_number
    for i, move in enumerate(pv):
        if board.turn == chess.WHITE:
            move_list_san.append(f"{move_num}. {board.san(move)}")
        else:
            move_list_san.append(board.san(move))
            move_num += 1
        board.push(move)

    final_sequence = " ".join(move_list_san)
    
    # Calculate mate in how many moves for White
    white_moves = (len(pv) + 1) // 2
    
    print(f"White has a forced checkmate in {white_moves} moves.")
    print("The winning sequence is:")
    print(final_sequence)
    
    # Close the engine
    await engine.quit()

if __name__ == "__main__":
    # In Python 3.7+
    asyncio.run(find_mate_in_n())