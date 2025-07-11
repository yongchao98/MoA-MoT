def solve_puzzle():
    """
    Solves the Fun Facts From The Zoo puzzle.
    1. Defines the answers to the riddles.
    2. Counts the occurrences of the letter 'A' in each answer.
    3. Prints the summation equation.
    4. Prints the final interpretation of the result.
    """
    answers = [
        "WHALE", "OUGHT", "REINDEER", "KOALA", "OCTOPUS", 
        "VIPER", "EEL", "RAT", "TUNA", "ANTLIA", "MANATEE", 
        "AARDVARK", "MOOSE", "NEWT", "ALPACA", "BEAR", 
        "EAGLE", "ANTEATER", "MACAQUE"
    ]

    # Step 2: Count the 'A's in each answer
    a_counts = [s.lower().count('a') for s in answers]

    # Step 3: Format and print the equation
    equation_str = " + ".join(map(str, a_counts))
    total_sum = sum(a_counts)
    
    print("The final equation is the sum of the number of 'A's in each answer word:")
    print(f"{equation_str} = {total_sum}")
    
    # Step 4: The final answer is derived from the sum
    # The sum is 21. A common three-word phrase related to 21 is "TWENTY ONE GUNS".
    
solve_puzzle()
