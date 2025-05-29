def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, occupied_positions):
    possible_moves = []
    directions = [(2, 1), (2, -1), (-2, 1), (-2, -1),
                  (1, 2), (1, -2), (-1, 2), (-1, -2)]
    for dx, dy in directions:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_position = new_x + new_y
        if 'A' <= new_x <= 'C' and '1' <= new_y <= '3' and new_position not in occupied_positions:
            possible_moves.append(new_position)
    return possible_moves

def solve_knight_swap():
    white_positions = ['A3', 'A1']
    black_positions = ['A2', 'B1']
    empty_positions = ['B3', 'C3', 'C2', 'C1']

    target_white_positions = ['A2', 'B1']
    target_black_positions = ['A3', 'A1']

    move_sequence = []

    # Move 1: w A3 to C2
    if is_valid_move('A3', 'C2'):
        white_positions.remove('A3')
        white_positions.append('C2')
        empty_positions.remove('C2')
        empty_positions.append('A3')
        move_sequence.append("w,A3,C2")

    # Move 2: B A2 to C3
    if is_valid_move('A2', 'C3'):
        black_positions.remove('A2')
        black_positions.append('C3')
        empty_positions.remove('C3')
        empty_positions.append('A2')
        move_sequence.append("B,A2,C3")

    # Move 3: w A1 to B3
    if is_valid_move('A1', 'B3'):
        white_positions.remove('A1')
        white_positions.append('B3')
        empty_positions.remove('B3')
        empty_positions.append('A1')
        move_sequence.append("w,A1,B3")

    # Move 4: B B1 to A1
    if is_valid_move('B1', 'A1'):
        black_positions.remove('B1')
        black_positions.append('A1')
        empty_positions.remove('A1')
        empty_positions.append('B1')
        move_sequence.append("B,B1,A1")

    # Move 5: w C2 to A2
    if is_valid_move('C2', 'A2'):
        white_positions.remove('C2')
        white_positions.append('A2')
        empty_positions.remove('A2')
        empty_positions.append('C2')
        move_sequence.append("w,C2,A2")

    # Move 6: B C3 to A3
    if is_valid_move('C3', 'A3'):
        black_positions.remove('C3')
        black_positions.append('A3')
        empty_positions.remove('A3')
        empty_positions.append('C3')
        move_sequence.append("B,C3,A3")

    if set(white_positions) == set(target_white_positions) and set(black_positions) == set(target_black_positions):
        return move_sequence
    else:
        return "No"

print(solve_knight_swap())