def is_on_diagonal(x, y):
    # Check if the point (x, y) is on the diagonal from (20, 1) to (1, 20)
    return x + y == 21

# Check if (14, 15) is on the diagonal
result = is_on_diagonal(14, 15)
print(result)