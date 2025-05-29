input_grid = [0, 2, 9, 6, 0, 4, 7, 0, 0, 8, 8, 0, 0, 0, 4, 2, 0, 2, 0, 7, 9, 0, 6, 6, 4, 0, 9, 0]
output_grid = [5 if x != 0 else 0 for x in input_grid]
print(output_grid)