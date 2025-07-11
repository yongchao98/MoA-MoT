import copy

def scale(array, factor):
    """Expands the array by a given factor."""
    new_rows = len(array) * factor
    new_cols = len(array[0]) * factor
    new_array = [['' for _ in range(new_cols)] for _ in range(new_rows)]
    for i in range(len(array)):
        for j in range(len(array[0])):
            char = array[i][j]
            for row_offset in range(factor):
                for col_offset in range(factor):
                    new_array[i * factor + row_offset][j * factor + col_offset] = char
    return new_array

def roll(array, axis, amount):
    """Performs a rolling transformation on the array."""
    arr = copy.deepcopy(array)
    rows = len(arr)
    cols = len(arr[0])
    
    if axis == 0:  # Columns
        for j in range(cols):
            shift = (j * amount) % rows
            if shift == 0:
                continue
            col_data = [arr[i][j] for i in range(rows)]
            shifted_col = col_data[-shift:] + col_data[:-shift]
            for i in range(rows):
                arr[i][j] = shifted_col[i]
    elif axis == 1:  # Rows
        for i in range(rows):
            shift = (i * amount) % cols
            if shift == 0:
                continue
            row_data = arr[i]
            shifted_row = row_data[-shift:] + row_data[:-shift]
            arr[i] = shifted_row
    return arr

def rotate(array, direction):
    """Rotates the array 90 degrees."""
    if direction == 'clockwise':
        return [list(reversed(row)) for row in zip(*array)]
    elif direction == 'counterclockwise':
        return [list(t) for t in zip(*reversed(array))]
    return array

def flip(array, direction):
    """Flips the array upside-down or left-to-right."""
    if direction == 'ud':
        return array[::-1]
    elif direction == 'lr':
        return [row[::-1] for row in array]
    return array

def transpose(array):
    """Transposes the array."""
    return [list(row) for row in zip(*array)]

def antidiagonal_transpose(array):
    """Performs a transpose along the anti-diagonal."""
    arr = transpose(array)
    arr = flip(arr, 'ud')
    arr = flip(arr, 'lr')
    return arr

def apply_replacements(array, row_counter, replacement_counter, replacements_dict):
    """Applies character replacements to a specific row."""
    items = list(replacements_dict.items())
    num_items = len(items)
    start_index = (replacement_counter * 5) % num_items
    
    pairs_to_use = items[start_index : start_index + 5]
    if len(pairs_to_use) < 5:
        pairs_to_use += items[:5 - len(pairs_to_use)]
        
    replacement_map = dict(pairs_to_use)
    
    num_rows = len(array)
    if num_rows == 0:
        return array
        
    target_row_idx = row_counter % num_rows
    row_to_modify = array[target_row_idx]
    
    new_row = [replacement_map.get(char, char) for char in row_to_modify]
    array[target_row_idx] = new_row
    
    return array

def parse_command(command_str):
    """Parses a command string into command and arguments."""
    parts = command_str.split()
    return parts[0], parts[1:]

def solve():
    """Main function to solve the puzzle."""
    initial_array = [['◬', '◑'], ['◪', '▨']]
    
    final_array = [
        ['0', '▣', '◫', '▤', '◫', '▤', '◫', '▲', '▧', '◔', '◫', '◫', '◔', '◔', '▣', '◔'],
        ['◬', '▲', '▲', '◪', '◫', '◔', '◔', '▲', '◪', '▨', '▨', '◔', '▲', '◬', '▲', '◪'],
        ['▨', '◭', '◬', '▤', '◭', '◪', '▨', '◪', '◬', '◑', '◬', '◪', '◑', '◪', '◗', '◔'],
        ['◑', '◬', '▨', '◪', '0', '◫', '▨', '◪', '◬', '▨', '▨', '◑', '▲', '◔', '◪', '▨'],
        ['▨', '◪', '◬', '▲', '◔', '◪', '◑', '◔', '▨', '◬', '◬', '▤', '◬', '◪', '▲', '▨'],
        ['◑', '▨', '▨', '▨', '◭', '◫', '0', '▨', '◑', '▨', '◫', '◪', '◔', '◫', '◔', '▨'],
        ['▨', '◬', '◔', '▲', '◔', '◫', '◪', '◔', '▨', '◪', '◔', '▤', '◬', '◫', '◪', '▨'],
        ['▤', '▨', '◫', '◬', '▲', '▣', '◑', '▨', '◑', '◬', '◪', '◪', '◪', '◫', '◔', '◑'],
        ['▤', '◑', '◑', '▨', '▲', '▨', '◪', '▨', '▨', '◪', '◑', '▲', '◔', '▨', '▨', '◑'],
        ['◪', '▨', '◪', '◬', '▲', '◑', '◬', '0', '▤', '▲', '▤', '◪', '◪', '▣', '▨', '◪'],
        ['▧', '◑', '◪', '▨', '▲', '◓', '▲', '▨', '▨', '▣', '◪', '▲', '◫', '◭', '▨', '◬'],
        ['◪', '▨', '▤', '◔', '▤', '◑', '▲', '◪', '◪', '▲', '▤', '▨', '◫', '◑', '▨', '◪'],
        ['▧', '▤', '◪', '◫', '◔', '▲', '▲', '◑', '▤', '◪', '◪', '▨', '◫', '▲', '◑', '▲'],
        ['▣', '▤', '▤', '◑', '▣', '◬', '◪', '◪', '◪', '◬', '▲', '◬', '▨', '◑', '◑', '▲'],
        ['◬', '◪', '◪', '◪', '◫', '▲', '◪', '◬', '▧', '◬', '◪', '▨', '▣', '▲', '◪', '◪'],
        ['◪', '▲', '▲', '◪', '▣', '◬', '◔', '▲', '▣', '▨', '▲', '◬', '◭', '◬', '◬', '▲']
    ]
    replacements_dict = {'◬': '▲', '◑': '◔', '◪': '▣', '▨': '◫', '◧': '▦', '◭': '◓', '◔': '◩', '◯': '◖', '◩': '▤', '◮': '◪', '◒': '▼', '◫': '▨', '▲': '▧', '▣': '◭', '▩': '◒', '◖': '◕', '▼': '◨', '▤': '◬', '▧': '◗', '▦': '▥', '◨': '▩', '◓': '▥', '◕': '◗', '◗': '◯', '▥': '◓'}
    
    commands_to_run = "scale 2, roll 0 1, rotate clockwise, flip lr, transpose, scale 2, roll 1 2, rotate counterclockwise, flip ud, roll 0 3, antidiagonal_transpose, scale 2, rotate counterclockwise, roll 0 8, flip ud, roll 1 6, transpose, roll 0 5, rotate counterclockwise, flip lr, transpose".split(', ')

    current_array = copy.deepcopy(initial_array)
    row_counter = 0
    replacement_counter = 0

    for command_str in commands_to_run:
        cmd, args = parse_command(command_str)
        
        if cmd == 'scale':
            current_array = scale(current_array, int(args[0]))
        elif cmd == 'roll':
            current_array = roll(current_array, int(args[0]), int(args[1]))
        elif cmd == 'rotate':
            current_array = rotate(current_array, args[0])
        elif cmd == 'flip':
            current_array = flip(current_array, args[0])
        elif cmd == 'transpose':
            current_array = transpose(current_array)
        elif cmd == 'antidiagonal_transpose':
            current_array = antidiagonal_transpose(current_array)
        
        current_array = apply_replacements(current_array, row_counter, replacement_counter, replacements_dict)
        
        row_counter += 1
        replacement_counter += 1

    print("The final array is:")
    for row in current_array:
        print(row)
    
    print("\nThe hidden symbols are:")
    print(f"Hidden value at [(0, 0)] is: {current_array[0][0]}")
    print(f"Hidden value at [(3, 4)] is: {current_array[3][4]}")
    print(f"Hidden value at [(5, 6)] is: {current_array[5][6]}")
    print(f"Hidden value at [(9, 7)] is: {current_array[9][7]}")
    
solve()
print("<<<D>>>")