# Testing different combinations to reach the target value of 24
def test_combinations():
    # Combination 1: (3 * 8) - (7 / 7)
    result1 = (3 * 8) - (7 / 7)
    
    # Combination 2: (3 * 7) + (7 - 8)
    result2 = (3 * 7) + (7 - 8)
    
    return result1, result2

# Print the results of the combinations
print(test_combinations())