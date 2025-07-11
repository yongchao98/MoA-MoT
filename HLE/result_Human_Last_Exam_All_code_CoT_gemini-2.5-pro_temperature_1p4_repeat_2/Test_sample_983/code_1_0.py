def find_next_number(n):
    """
    Calculates the next number in the sequence based on the rule:
    Next = Current + (sum of all digits) * (sum of middle two digits).
    """
    if n < 1000 or n > 9999:
        print("The number must be a 4-digit number.")
        return

    s_n = str(n)
    d1 = int(s_n[0])
    d2 = int(s_n[1])
    d3 = int(s_n[2])
    d4 = int(s_n[3])

    sum_all_digits = d1 + d2 + d3 + d4
    sum_middle_digits = d2 + d3
    
    difference = sum_all_digits * sum_middle_digits
    next_num = n + difference
    
    print(f"The next number is found by the calculation:")
    print(f"{n} + ({d1} + {d2} + {d3} + {d4}) * ({d2} + {d3}) = {n} + {sum_all_digits} * {sum_middle_digits} = {n} + {difference} = {next_num}")
    
    # Return the final number for the closing tag.
    return next_num

# The last number in the provided sequence
last_number = 2352
final_answer = find_next_number(last_number)
