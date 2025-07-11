def solve_riddle():
    """
    This function solves the number riddle based on a step-by-step logical deduction
    and prints the final sequence.
    """
    
    # Initialize a list of 5 positions, using 0 as a placeholder.
    sequence = [0, 0, 0, 0, 0]
    
    print("Let's solve the riddle step by step:")
    
    # Clue: "...the fifth, who never had a thing and lastly follows."
    # This places the number 5 in the last position (index 4).
    sequence[4] = 5
    print(f"1. The fifth (5) 'lastly follows', so it is at the end: {sequence}")
    
    # Clue: "Number three is the strongest and it protects the last;"
    # This places the number 3 just before the last one, in the fourth position (index 3).
    sequence[3] = 3
    print(f"2. Number 3 'protects the last', so it is in the fourth position: {sequence}")
    
    # Clue: "Number 1 is the first to have something..."
    # This places the number 1 in the first position (index 0).
    sequence[0] = 1
    print(f"3. Number 1 'is the first', so it goes in the first position: {sequence}")
    
    # Clue: "Number 4 ... always goes behind [Number 2]."
    # The remaining numbers are 2 and 4. The remaining positions are the second and third.
    # For 4 to be 'behind' 2, 2 must come first.
    sequence[1] = 2
    sequence[2] = 4
    print(f"4. Number 4 must go 'behind' 2. With the remaining slots, 2 must be second and 4 third: {sequence}")
    
    # Final Result
    print("\nThe final sequence is formed by the numbers in their determined positions.")
    final_equation = ", ".join(map(str, sequence))
    print(f"The final order is: {final_equation}")

solve_riddle()