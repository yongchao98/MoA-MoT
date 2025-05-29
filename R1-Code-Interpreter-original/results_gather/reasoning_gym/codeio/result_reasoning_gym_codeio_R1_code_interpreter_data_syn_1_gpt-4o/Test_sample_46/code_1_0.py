import itertools

def find_largest_group(matrix):
    def dfs(x, y, part_number):
        stack = [(x, y)]
        count = 0
        while stack:
            cx, cy = stack.pop()
            if 0 <= cx < len(matrix) and 0 <= cy < len(matrix[0]) and matrix[cx][cy] == part_number:
                count += 1
                matrix[cx][cy] = -1  # Mark as visited
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    stack.append((cx + dx, cy + dy))
        return count

    max_group_size = 0
    part_number = 0

    for i, j in itertools.product(range(len(matrix)), range(len(matrix[0]))):
        if matrix[i][j] > 0:
            group_size = dfs(i, j, matrix[i][j])
            if group_size > max_group_size:
                max_group_size = group_size
                part_number = matrix[i][j]

    return [max_group_size, part_number]

# Test with a matrix filled with zeros
matrix = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

result = find_largest_group(matrix)
print(result)