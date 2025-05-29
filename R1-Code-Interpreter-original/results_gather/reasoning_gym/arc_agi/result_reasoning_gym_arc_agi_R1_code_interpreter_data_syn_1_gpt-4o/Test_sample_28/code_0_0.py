def extract_non_six_numbers(input_grid):
    non_six_numbers = []
    for row in input_grid:
        for num in row:
            if num != 6:
                non_six_numbers.append(num)
    return non_six_numbers

def form_output_grid(non_six_numbers):
    # Assuming the order is based on the examples, we fill the 3x3 grid
    return [
        [non_six_numbers[0], non_six_numbers[1], non_six_numbers[2]],
        [non_six_numbers[3], non_six_numbers[4], non_six_numbers[5]],
        [non_six_numbers[6], non_six_numbers[7], non_six_numbers[8]]
    ]

# Test input grid
input_grid = [
    [6, 6, 0, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 4, 6, 6, 6],
    [3, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 1, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 0, 6],
    [6, 5, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 3, 6],
    [6, 6, 6, 6, 6, 2, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 9, 6, 6]
]

non_six_numbers = extract_non_six_numbers(input_grid)
output_grid = form_output_grid(non_six_numbers)

print(output_grid)