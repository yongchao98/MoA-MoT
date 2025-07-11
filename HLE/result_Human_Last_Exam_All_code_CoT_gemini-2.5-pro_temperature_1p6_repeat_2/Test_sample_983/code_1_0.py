def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    while n:
        s += n % 10
        n //= 10
    return s

def find_next_number(sequence):
    """
    Finds the next number in the sequence.
    The problem is very complex, and a single rule does not seem to fit the entire sequence.
    The most evident pattern is at the beginning (differences of 6, 12, 18), which suggests
    the difference increases by 6 each time (an arithmetic progression of differences).
    We will apply this logic to the last known step.
    """
    last_number = sequence[-1]
    second_last_number = sequence[-2]
    
    # Calculate the last difference
    last_difference = last_number - second_last_number
    
    # Assume the difference of differences is 6
    next_difference = last_difference + 6
    
    # Calculate the next number in the sequence
    next_number = last_number + next_difference
    
    # Print the equation
    print(f"The last number is {last_number}.")
    print(f"The previous number is {second_last_number}.")
    print(f"The last difference is {last_number} - {second_last_number} = {last_difference}.")
    print(f"Assuming the difference of differences is 6, the next difference is {last_difference} + 6 = {next_difference}.")
    print(f"So, the next number in the sequence is {last_number} + {next_difference} = {next_number}.")
    return next_number

# The given sequence
sequence = [1014, 1020, 1032, 1050, 1104, 1224, 1230, 1290, 1410, 1440, 1452, 1482, 1500, 1512, 1560, 1584, 1602, 1734, 1812, 1890, 1902, 1932, 1974, 2004, 2040, 2064, 2070, 2130, 2142, 2274, 2340, 2352]

result = find_next_number(sequence)