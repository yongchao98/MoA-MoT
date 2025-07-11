def solve():
    """
    This script identifies the person from the riddle and then performs
    a numerological calculation on their last name as a coding example.
    """
    # The last name of the man described in the riddle.
    last_name = "Rezanov"

    # Convert the name to uppercase for consistent letter values.
    name_upper = last_name.upper()

    # List to hold the numerical value of each letter.
    values = []
    # List to hold the string parts of the equation, e.g., "R(18)".
    equation_parts = []

    # Iterate through each character in the name.
    for char in name_upper:
        if 'A' <= char <= 'Z':
            # Calculate the value (A=1, B=2, ... Z=26) using ASCII codes.
            # ord('A') is 65. So for 'A', it's 65 - 65 + 1 = 1.
            # For 'B', it's 66 - 65 + 1 = 2.
            value = ord(char) - ord('A') + 1
            values.append(value)
            equation_parts.append(f"{char}({value})")

    # Calculate the total sum.
    total = sum(values)

    # Join the parts of the equation with a " + ".
    equation_string = " + ".join(equation_parts)

    print(f"The last name is {last_name}.")
    print("Here is a calculation based on the letters in the name:")
    # Print the final equation with each number and the sum.
    print(f"{equation_string} = {total}")

solve()