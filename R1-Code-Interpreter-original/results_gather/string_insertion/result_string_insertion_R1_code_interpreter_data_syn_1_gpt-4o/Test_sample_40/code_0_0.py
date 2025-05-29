def apply_t_operation(s):
    i = 0
    result = list(s)
    while i < len(result) - 3:
        if result[i:i+4] == ['A', 'B', 'C', 'D']:
            result.insert(i+4, 'A')
            i += 5  # Move past the inserted character
        elif result[i:i+4] == ['B', 'C', 'D', 'E']:
            result.insert(i+4, 'B')
            i += 5
        elif result[i:i+4] == ['C', 'D', 'E', 'A']:
            result.insert(i+4, 'C')
            i += 5
        elif result[i:i+4] == ['D', 'E', 'A', 'B']:
            result.insert(i+4, 'D')
            i += 5
        elif result[i:i+4] == ['E', 'A', 'B', 'C']:
            result.insert(i+4, 'E')
            i += 5
        else:
            i += 1
    return ''.join(result)

# Given string
s = "CEABCCBCDEDCDEAAD"
# Apply the T operation
result = apply_t_operation(s)
print(result)