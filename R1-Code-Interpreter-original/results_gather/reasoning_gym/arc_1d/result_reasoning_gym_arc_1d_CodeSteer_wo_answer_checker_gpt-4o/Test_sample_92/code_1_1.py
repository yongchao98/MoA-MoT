from collections import deque

def transform_grid(input_grid):
    queue = deque()
    output_grid = []

    for num in input_grid:
        if num != 0:
            queue.append(num)
            output_grid.append(num)
        else:
            if queue:
                output_grid.append(queue[0])
            else:
                output_grid.append(0)

    output_grid.append(0)  # Append a zero at the end
    return output_grid

# Test input
input_grid = [8, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 5, 2, 0, 0, 0, 3, 0, 8, 5, 0, 0, 0, 0, 0, 2, 8]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")