def transform_grid(input_grid):
    output_grid = []
    in_sequence = False

    for num in input_grid:
        if num == 2:
            if not in_sequence:
                in_sequence = True
            output_grid.append(8)
        elif num == 0:
            in_sequence = False
            output_grid.append(0)
        else:
            output_grid.append(num)

    return output_grid

# Test input
input_grid = [0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")