# Trying different combinations to reach 654
def find_combination():
    numbers = [3, 3, 61, 36, 75]
    
    # Try multiplying different pairs and adding/subtracting the rest
    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                # Multiply two numbers
                product = numbers[i] * numbers[j]
                
                # Try adding the remaining numbers
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                
                # Check all combinations of adding/subtracting remaining numbers
                for a in remaining:
                    for b in remaining:
                        if a != b:
                            if product + a + b == 654:
                                return f"{numbers[i]} * {numbers[j]} + {a} + {b}"
                            if product + a - b == 654:
                                return f"{numbers[i]} * {numbers[j]} + {a} - {b}"
                            if product - a + b == 654:
                                return f"{numbers[i]} * {numbers[j]} - {a} + {b}"
                            if product - a - b == 654:
                                return f"{numbers[i]} * {numbers[j]} - {a} - {b}"
    
    return "No combination found"

print(find_combination())