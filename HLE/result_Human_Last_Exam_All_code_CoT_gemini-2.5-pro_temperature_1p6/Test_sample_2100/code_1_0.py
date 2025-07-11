def solve():
    """
    This function calculates a numerical value for the name "Rezanov"
    by summing the alphabetical positions of its letters (A=1, B=2, ...).
    It then prints the full equation.
    """
    name = "Rezanov"
    name_upper = name.upper()
    
    numbers = []
    # Loop through each character in the name
    for char in name_upper:
        # Check if the character is an alphabet letter
        if 'A' <= char <= 'Z':
            # Calculate the position in the alphabet (A=1, B=2, etc.)
            # ord(char) gets the ASCII value of the character.
            # ord('A') is the ASCII value of 'A'.
            # The difference + 1 gives the alphabetical position.
            value = ord(char) - ord('A') + 1
            numbers.append(value)
            
    # Calculate the sum of all numbers
    total = sum(numbers)
    
    # Create the string representation of the equation
    equation_parts = [str(n) for n in numbers]
    equation = " + ".join(equation_parts)
    
    # Print the final result in the format: 18 + 5 + ... = 101
    print(f"{equation} = {total}")

solve()