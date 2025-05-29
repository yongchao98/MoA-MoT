# Evaluate different combinations to find the one that results in 24
def evaluate_combinations():
    results = []
    # Combination 1: (7 + 3) * (4 + 8)
    result1 = (7 + 3) * (4 + 8)
    results.append(result1)
    
    # Combination 2: (7 * 3) + (4 * 8)
    result2 = (7 * 3) + (4 * 8)
    results.append(result2)
    
    # Combination 3: 7 + (3 * 4) + 8
    result3 = 7 + (3 * 4) + 8
    results.append(result3)
    
    # Combination 4: 7 * (3 + 4) - 8
    result4 = 7 * (3 + 4) - 8
    results.append(result4)
    
    return results

print(evaluate_combinations())