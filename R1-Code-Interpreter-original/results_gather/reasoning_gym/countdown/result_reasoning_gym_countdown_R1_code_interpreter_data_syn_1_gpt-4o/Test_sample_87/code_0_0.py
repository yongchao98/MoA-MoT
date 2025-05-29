# Trying different combinations to reach 654
def find_combination():
    numbers = [3, 3, 61, 36, 75]
    
    # Try multiplying 75 by 3
    result1 = 75 * 3
    
    # Check if we can reach 654 by adding or subtracting the other numbers
    for num in numbers:
        if num != 75 and num != 3:
            result2 = result1 * num
            if result2 == 654:
                return f"75 * 3 * {num}"
    
    # Try other combinations
    result3 = 75 * 36
    if result3 == 654:
        return "75 * 36"
    
    result4 = 75 * 3 * 3
    if result4 == 654:
        return "75 * 3 * 3"
    
    # Try adding or subtracting
    result5 = 75 * 3 + 36
    if result5 == 654:
        return "75 * 3 + 36"
    
    result6 = 75 * 3 + 61
    if result6 == 654:
        return "75 * 3 + 61"
    
    return "No combination found"

print(find_combination())