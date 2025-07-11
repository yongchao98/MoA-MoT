def solve_medallion_dating():
    """
    This function determines and prints the century of the Ténès medallion.
    """
    # Based on archaeological dating, the medallion is from the 4th century A.D.
    century_number = 4

    # The instruction asks to output each number in the final equation.
    # As there is no complex equation, I will output the key number used for the answer.
    print(f"The number representing the century is: {century_number}")

    # Determine the correct ordinal suffix for the century number.
    if 11 <= (century_number % 100) <= 13:
        suffix = 'th'
    else:
        last_digit = century_number % 10
        if last_digit == 1:
            suffix = 'st'
        elif last_digit == 2:
            suffix = 'nd'
        elif last_digit == 3:
            suffix = 'rd'
        else:
            suffix = 'th'
    
    # Define the era.
    era = "A.D."

    # Construct the final answer string. While the example format is "9 A.D.",
    # "4th century A.D." is the correct and unambiguous way to state a century.
    final_answer = f"{century_number}{suffix} century {era}"

    # Print the final answer.
    print(f"Archaeologists ascribe this medallion to the {final_answer}")

solve_medallion_dating()