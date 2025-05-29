from itertools import permutations

def compute_all_expressions(numbers):
    if len(numbers) == 1:
        return numbers  # Already in the form (value, expression)
    
    results = []
    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                # Pick two numbers to operate on
                num1, expr1 = numbers[i]
                num2, expr2 = numbers[j]
                # Remaining numbers
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                
                # Try all operations
                results.extend(compute_all_expressions([(num1 + num2, f"({expr1} + {expr2})")] + remaining))
                results.extend(compute_all_expressions([(num1 - num2, f"({expr1} - {expr2})")] + remaining))
                results.extend(compute_all_expressions([(num1 * num2, f"({expr1} * {expr2})")] + remaining))
                if num2 != 0:
                    results.extend(compute_all_expressions([(num1 / num2, f"({expr1} / {expr2})")] + remaining))
    
    return results

def find_expression_for_24(numbers):
    for perm in permutations(numbers):
        # Convert numbers to tuples of (value, expression)
        initial = [(num, str(num)) for num in perm]
        results = compute_all_expressions(initial)
        for result, expression in results:
            if abs(result - 24) < 1e-6:  # Use a tolerance for floating-point comparison
                return expression
    return None

numbers = [1, 5, 10, 12]
expression = find_expression_for_24(numbers)
if expression:
    print(f"<<<{expression} = 24>>>")
else:
    print("No valid expression found.")