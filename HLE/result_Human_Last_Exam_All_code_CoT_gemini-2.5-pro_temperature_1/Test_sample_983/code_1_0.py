def get_digits(n):
    """Returns a list of digits of a number."""
    return [int(d) for d in str(n)]

def solve_sequence(n):
    """
    Calculates the next number in the sequence based on the rule:
    next = current + (sum of digits) * (largest digit)
    """
    if n < 0:
        n = abs(n)
        
    digits = get_digits(n)
    
    # Handle the case of single-digit numbers
    if not digits:
        return 0
        
    sum_of_digits = sum(digits)
    largest_digit = max(digits)
    
    increment = sum_of_digits * largest_digit
    next_number = n + increment
    
    print(f"The last number in the sequence is {n}.")
    print(f"The sum of its digits is {' + '.join(map(str, digits))} = {sum_of_digits}.")
    print(f"The largest digit is {largest_digit}.")
    print(f"The increment is calculated as {sum_of_digits} * {largest_digit} = {increment}.")
    print(f"The next number in the sequence is {n} + {increment} = {next_number}.")
    
    return next_number

# The last number in the provided sequence.
last_number = 2352
final_answer = solve_sequence(last_number)
# The final answer is printed within the function.