def find_cantor_pair(target):
    for a in range(1000):  # Arbitrary range, can be adjusted
        for b in range(1000):  # Arbitrary range, can be adjusted
            result = 0.5 * (a + b) * (a + b + 1) + b
            if result == target:
                return a, b
    return None

# Given output
target_value = 127700.0

# Find the pair (a, b)
pair = find_cantor_pair(target_value)
print(pair)