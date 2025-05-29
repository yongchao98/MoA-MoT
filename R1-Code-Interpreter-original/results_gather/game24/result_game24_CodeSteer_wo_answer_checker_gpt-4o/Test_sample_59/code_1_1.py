from itertools import permutations

def compute_all_expressions(numbers):
    if len(numbers) == 1:
        return [numbers[0]]
    
    results = []
    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                # Pick two numbers to operate on
                num1 = numbers[i]
                num2 = numbers[j]
                # Remaining numbers
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                
                # Try all operations
                results.extend(compute_all_expressions([num1 + num2] + remaining))
                results.extend(compute_all_expressions([num1 - num2] + remaining))
                results.extend(compute_all_expressions([num1 * num2] + remaining))
                if num2 != 0:
                    results.extend(compute_all_expressions([num1 / num2] + remaining))
    
    return results

def find_expression_for_24(numbers):
    for perm in permutations(numbers):
        results = compute_all_expressions(list(perm))
        if 24 in results:
            return True
    return False

numbers = [1, 5, 10, 12]
if find_expression_for_24(numbers):
    print("A valid expression that evaluates to 24 was found.")
else:
    print("No valid expression found.")

# Note: This code checks if a valid expression exists but does not output the exact expression.
# To output the exact expression, additional logic is needed to track the operations and format the expression.