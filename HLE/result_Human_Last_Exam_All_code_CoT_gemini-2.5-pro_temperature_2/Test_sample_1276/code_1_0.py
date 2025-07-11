def solve_zoo_puzzle():
    """
    Solves the 'Fun Facts From The Zoo' puzzle by calculating letter values.
    """
    animals = [
        "HUMAN", "OTTER", "REINDEER", "KOALA", "OCTOPUS", "VIPER", "EEL", "RAT",
        "TUNA", "ANT", "NAUTILUS", "RHINO", "LEMMING", "NEWT", "ALPACA", "TOAD",
        "EAGLE", "ANTEATER", "APE"
    ]
    numbers = [5, 5, 8, 5, 7, 5, 3, 3, 4, 6, 7, 8, 5, 4, 6, 4, 5, 8, 7]
    
    final_phrase = []
    
    print("Executing calculation: (sum of letter values of animal name) + number = result -> letter\n")

    for i, animal in enumerate(animals):
        # Calculate sum of letter values (A=1, B=2, ...)
        animal_sum = sum(ord(char.upper()) - ord('A') + 1 for char in animal)
        
        # Get the number from the puzzle
        num_from_puzzle = numbers[i]
        
        # Perform the final calculation
        total = animal_sum + num_from_puzzle
        
        # Convert back to a letter (1-26 -> A-Z)
        # (total - 1) % 26 maps the result to 0-25
        # + 1 maps it back to 1-26
        final_letter_num = (total - 1) % 26 + 1
        final_letter = chr(final_letter_num + ord('A') - 1)
        final_phrase.append(final_letter)
        
        print(f"Riddle {i+1}: {animal} ({animal_sum}) + {num_from_puzzle} = {total} => {final_letter}")

    print("\nFinal decoded phrase:")
    print("".join(final_phrase))

solve_zoo_puzzle()
<<<JUST BE YOURSELF>>>