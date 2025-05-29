# Mapping from coordinates to chessboard notation
def coord_to_chess_notation(coord):
    col_map = {0: 'A', 1: 'B', 2: 'C', 3: 'D'}
    row_map = {0: '1', 1: '2', 2: '3'}
    x, y = coord
    return f"{col_map[y]}{row_map[x]}"

# Convert the output to the required format
def convert_output(moves):
    formatted_moves = []
    for move in moves:
        # Extract color, from_coord, and to_coord using string manipulation
        color, from_coord, to_coord = move[0], move[2:7], move[8:13]
        from_chess = coord_to_chess_notation(eval(from_coord))
        to_chess = coord_to_chess_notation(eval(to_coord))
        formatted_moves.append(f"{color},{from_chess},{to_chess}")
    return formatted_moves

# Original output from the code
original_moves = ['w,(1, 2),(2, 0)', 'B,(2, 1),(1, 3)', 'w,(2, 0),(1, 2)', 'B,(1, 3),(1, 3)', 
                  'w,(1, 2),(1, 2)', 'B,(1, 3),(0, 1)', 'w,(1, 2),(0, 0)', 'B,(0, 1),(0, 1)', 
                  'w,(0, 0),(0, 0)', 'B,(0, 1),(2, 0)', 'w,(0, 0),(2, 1)', 'B,(2, 0),(2, 0)', 
                  'w,(2, 1),(2, 1)', 'B,(2, 0),(1, 2)']

# Convert and print the formatted moves
formatted_moves = convert_output(original_moves)
print(formatted_moves)