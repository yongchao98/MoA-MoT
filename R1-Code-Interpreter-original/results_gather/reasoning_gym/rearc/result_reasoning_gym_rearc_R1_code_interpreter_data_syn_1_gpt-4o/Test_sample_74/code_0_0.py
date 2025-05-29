def find_distinct_numbers(grid):
    seen = set()
    result = []
    for row in grid:
        for num in row:
            if num not in seen:
                seen.add(num)
                result.append(num)
    return result

input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 5, 0, 0, 0, 0, 0],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = find_distinct_numbers(input_grid)
for num in output_grid:
    print(num)