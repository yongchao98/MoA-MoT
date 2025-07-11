import copy

def solve_puzzle():
    """
    This function sets up the puzzle, runs the correct sequence of commands,
    and prints the final state and the hidden symbols.
    """
    initial_array = [
        ['◬', '◑'],
        ['◪', '▨']
    ]

    replacements_dict = {
        '◬': '▲', '◑': '◔', '◪': '▣', '▨': '◫', '◧': '▦', '◭': '◓', '◔': '◩', '◯': '◖',
        '◩': '▤', '◮': '◪', '◒': '▼', '◫': '▨', '▲': '▧', '▣': '◭', '▩': '◒', '◖': '◕',
        '▼': '◨', '▤': '◬', '▧': '◗', '▦': '▥', '◨': '▩', '◓': '▥', '◕': '◗', '◗': '◯',
        '▥': '◓'
    }

    # --- Transformation Functions ---
    def scale(array, factor):
        new_rows, new_cols = len(array) * factor, len(array[0]) * factor
        new_array = [['' for _ in range(new_cols)] for _ in range(new_rows)]
        for r_idx, row in enumerate(array):
            for c_idx, char in enumerate(row):
                for i in range(factor):
                    for j in range(factor):
                        new_array[r_idx * factor + i][c_idx * factor + j] = char
        return new_array

    def roll(array, axis, amount):
        arr = copy.deepcopy(array)
        num_rows, num_cols = len(arr), len(arr[0])
        if axis == 0:  # Columns shift downward
            new_arr = copy.deepcopy(arr)
            for j in range(num_cols):
                shift = j * amount
                for i in range(num_rows):
                    new_arr[i][j] = arr[(i - shift) % num_rows][j]
            return new_arr
        elif axis == 1:  # Rows shift right
            new_arr = copy.deepcopy(arr)
            for i in range(num_rows):
                shift = i * amount
                original_row = arr[i]
                new_arr[i] = original_row[-shift % num_cols:] + original_row[:-shift % num_cols]
            return new_arr
        return arr

    def rotate(array, direction):
        rows, cols = len(array), len(array[0])
        new_array = [['' for _ in range(rows)] for _ in range(cols)]
        if direction == 'clockwise':
            for r in range(rows):
                for c in range(cols):
                    new_array[c][rows - 1 - r] = array[r][c]
        elif direction == 'counterclockwise':
            for r in range(rows):
                for c in range(cols):
                    new_array[cols - 1 - c][r] = array[r][c]
        return new_array

    def flip(array, direction):
        if direction == 'ud':
            return array[::-1]
        elif direction == 'lr':
            return [row[::-1] for row in array]
        return array

    def transpose(array):
        rows, cols = len(array), len(array[0])
        new_array = [['' for _ in range(rows)] for _ in range(cols)]
        for r in range(rows):
            for c in range(cols):
                new_array[c][r] = array[r][c]
        return new_array

    def antidiagonal_transpose(array):
        rows, cols = len(array), len(array[0])
        new_array = [['' for _ in range(rows)] for _ in range(cols)]
        for r in range(rows):
            for c in range(cols):
                new_array[cols - 1 - c][rows - 1 - r] = array[r][c]
        return new_array

    # --- Main Processing Logic ---
    command_str = "scale 2, roll 0 1, rotate clockwise, flip lr, transpose, scale 2, roll 1 2, rotate counterclockwise, flip ud, roll 0 3, antidiagonal_transpose, scale 2, rotate counterclockwise, roll 0 8, flip ud, roll 1 5, transpose, roll 0 5, rotate counterclockwise, flip lr, transpose"
    commands = [cmd.strip() for cmd in command_str.split(',')]
    
    array = copy.deepcopy(initial_array)
    row_counter = 0
    replacement_offset = 0
    all_keys = list(replacements_dict.keys())

    for cmd_str in commands:
        parts = cmd_str.split()
        command = parts[0]

        if command == 'scale':
            array = scale(array, int(parts[1]))
        elif command == 'roll':
            array = roll(array, int(parts[1]), int(parts[2]))
        elif command == 'rotate':
            array = rotate(array, parts[1])
        elif command == 'flip':
            array = flip(array, parts[1])
        elif command == 'transpose':
            array = transpose(array)
        elif command == 'antidiagonal_transpose':
            array = antidiagonal_transpose(array)
        
        # Apply replacements
        num_rows = len(array)
        target_row_index = row_counter % num_rows
        keys_to_replace = [all_keys[i % len(all_keys)] for i in range(replacement_offset, replacement_offset + 5)]

        for j in range(len(array[target_row_index])):
            char = array[target_row_index][j]
            if char in keys_to_replace:
                array[target_row_index][j] = replacements_dict[char]

        row_counter += 1
        replacement_offset = (replacement_offset + 5) % len(all_keys)

    # --- Output Results ---
    print("Final calculated array:")
    for row in array:
        print(row)
    
    print("\nHidden Symbols:")
    hidden_coords = [(0, 0), (3, 4), (5, 6), (9, 7)]
    for r, c in hidden_coords:
        print(f"Hidden value at [({r}, {c})] is: {array[r][c]}")

solve_puzzle()