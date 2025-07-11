def get_digits(n):
    """Helper function to get the four digits of a number."""
    d3 = n // 1000
    d2 = (n // 100) % 10
    d1 = (n // 10) % 10
    d0 = n % 10
    return d3, d2, d1, d0

def calculate_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    next_num = current_num + 6 * ((d3*d1) + (d2*d0))
    """
    d3, d2, d1, d0 = get_digits(n)
    
    # Calculate the increment k based on the proposed rule
    k = (d3 * d1) + (d2 * d0)
    
    # The total increment is 6 * k
    increment = 6 * k
    
    # Calculate the next number
    next_num = n + increment
    
    print(f"The last number in the sequence is: {n}")
    print(f"Let the four digits of {n} be a, b, c, d, where a={d3}, b={d2}, c={d1}, d={d0}")
    print(f"The discovered rule for the increment is: 6 * (a*c + b*d)")
    print(f"Applying the rule: 6 * ({d3}*{d1} + {d2}*{d0}) = 6 * {k} = {increment}")
    print(f"The final equation is: {n} + {increment} = {next_num}")
    print(f"So, the next number in the sequence is: {next_num}")
    return next_num

# The last number of the given sequence
last_number = 2352

# Calculate and print the next number
final_answer = calculate_next_number(last_number)