def is_length_perfect_square(w: str):
    """
    This function implements a decision procedure to determine if the length
    of a string 'w' is a perfect square. It is guaranteed to halt for any input.
    
    Args:
        w: The input string.
    
    Returns:
        None. It prints the decision process and result.
    """
    n = len(w)
    print(f"Input string: \"{w}\"")
    print(f"Length of string is n = {n}")
    print(f"Starting check to see if n = {n} is a perfect square.")
    print("-" * 20)

    # The decider loop
    i = 0
    while True:
        square = i * i
        
        if square == n:
            # Found a perfect square match
            print(f"Checking i={i}: {i} * {i} = {square}")
            print("\nResult: Found a match!")
            print(f"The length {n} is a perfect square.")
            print(f"The final equation is: {i} * {i} == {n}")
            break
        elif square > n:
            # Overshot n, so it cannot be a perfect square
            print(f"Checking i={i}: {i} * {i} = {square}")
            print(f"\nResult: Stop! {square} is greater than {n}.")
            print(f"The length {n} is not a perfect square.")
            break
        else:
            # Keep searching
            print(f"Checking i={i}: {i} * {i} = {square} (no match, continue...)")
            i += 1
            
    print("-" * 20)
    print("The decision procedure has halted.")


# Main part of the script
# Let's test with a string whose length is a perfect square (e.g., 25)
input_string = "123456789-123456789-12345"
is_length_perfect_square(input_string)

print("\nBecause we can construct an algorithm that always halts and correctly decides membership,")
print("the language L = {w | |w| is a perfect square} is decidable.")
