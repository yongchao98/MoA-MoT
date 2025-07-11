def solve():
    """
    Calculates the next number in the sequence based on the names of
    chemical elements.
    The rule is the sum of the squares of the letters' 1-based
    alphabetical positions.
    """
    word = "Selenium"  # The 34th element, following the 33 numbers given.
    
    total = 0
    equation_parts = []
    
    # Calculate the value for the word
    for char in word:
        # Get 1-based position in the alphabet (a=1, b=2, ...)
        pos = ord(char.lower()) - ord('a') + 1
        square = pos * pos
        total += square
        equation_parts.append(f"{pos}^2")

    # Print the full equation as requested
    equation_str = " + ".join(equation_parts)
    print(f"The calculation for '{word}' is based on the sum of the squares of its letters' positions in the alphabet (A=1, Z=26).")
    print(f"Calculation: {equation_str} = {total}")
    print(f"\nThe next number in the sequence is: {total}")

solve()