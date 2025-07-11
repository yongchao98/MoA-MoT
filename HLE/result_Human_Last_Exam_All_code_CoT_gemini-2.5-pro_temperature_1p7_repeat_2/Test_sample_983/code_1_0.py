def find_next_number_in_sequence():
    """
    Calculates the next number in a sequence based on a derived rule.
    The rule is: next = current + 6 * (d1*d3 + d2*d4), where the current number
    is represented by the digits d1, d2, d3, d4.
    """
    sequence = [
        1014, 1020, 1032, 1050, 1104, 1224, 1230, 1290, 1410, 1440,
        1452, 1482, 1500, 1512, 1560, 1584, 1602, 1734, 1812, 1890,
        1902, 1932, 1974, 2004, 2040, 2064, 2070, 2130, 2142, 2274,
        2340, 2352
    ]
    
    last_number = sequence[-1]
    
    # Extract digits from the last number
    s_num = str(last_number)
    if len(s_num) != 4:
        print("The logic is designed for 4-digit numbers.")
        return

    d1 = int(s_num[0])
    d2 = int(s_num[1])
    d3 = int(s_num[2])
    d4 = int(s_num[3])
    
    # Apply the derived rule
    increment_factor = (d1 * d3) + (d2 * d4)
    increment = 6 * increment_factor
    next_number = last_number + increment
    
    # Print the final equation and the result
    print(f"{last_number} + 6 * ({d1} * {d3} + {d2} * {d4}) = {next_number}")
    print(f"The next number in the sequence is: {next_number}")

find_next_number_in_sequence()
<<<2448>>>