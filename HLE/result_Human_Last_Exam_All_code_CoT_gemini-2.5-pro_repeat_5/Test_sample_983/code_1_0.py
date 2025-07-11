def get_digits(n):
    """Extracts the four digits from a four-digit number."""
    return [int(d) for d in str(n).zfill(4)]

def calculate_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    a_{n+1} = a_n + 6 * k, where k = (d1 * d2) + d3.
    """
    if n < 1000 or n > 9999:
        # The rule is defined for 4-digit numbers.
        # This part is for safety, assuming the pattern continues.
        s = str(n)
        d1 = int(s[0])
        d2 = int(s[1])
        d3 = int(s[2])
    else:
      digits = get_digits(n)
      d1, d2, d3, _ = digits

    # Calculate k based on the rule
    k = (d1 * d2) + d3
    
    # Calculate the difference to be added
    diff = 6 * k
    
    # Calculate the next number
    next_num = n + diff
    
    print(f"The last number is {n}.")
    print(f"The digits are: d1={d1}, d2={d2}, d3={d3}.")
    print(f"The value of k is calculated as (d1 * d2) + d3 = ({d1} * {d2}) + {d3} = {k}.")
    print(f"The difference is 6 * k = 6 * {k} = {diff}.")
    print(f"The next number in the sequence is {n} + {diff} = {next_num}.")
    return next_num

# The last number in the given sequence
last_number = 2352

# Calculate and print the next number
next_number_in_sequence = calculate_next_number(last_number)