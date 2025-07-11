import math

def solve_letter_arrangement_puzzle():
    """
    Calculates the number of ways to arrange the letters L, N, S, W
    based on the given connection rule.
    """
    
    # Based on the analysis, only letters ending in a "Low" position can be
    # followed by another letter. These letters are {L, N, S}.
    num_arrangeable_letters = 3
    
    # The number of ways to arrange these letters is the factorial of the count.
    result = math.factorial(num_arrangeable_letters)
    
    # Create the equation string as requested, e.g., "3 * 2 * 1"
    equation_parts = list(range(num_arrangeable_letters, 0, -1))
    equation_str = " * ".join(map(str, equation_parts))
    
    print("Based on the connection rules, the letter 'W' must be the last letter in the sequence.")
    print("The other three letters {L, N, S} can be arranged in any order before 'W'.")
    print(f"The number of ways to arrange these {num_arrangeable_letters} letters is {num_arrangeable_letters}! (factorial).")
    print("The calculation is:")
    
    # Print the full equation including each number and the result.
    print(f"{equation_str} = {result}")

solve_letter_arrangement_puzzle()