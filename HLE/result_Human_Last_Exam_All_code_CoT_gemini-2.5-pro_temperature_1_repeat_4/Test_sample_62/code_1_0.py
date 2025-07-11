def solve_riddle():
    """
    Solves the riddle based on the clues from Elizabeth George's novel.
    """
    print("Decoding the riddle about Christian's preference...")
    print("-" * 40)

    # Clue 1: The subject "THEM" is types of coffee.
    print("Step 1: The clue 'Cafi...' is a mispronunciation of the Italian word for coffee.")
    first_part = "Caffè"
    print(f"The first part of the name is: {first_part}")

    # Clue 2: The preference for a "peculiar shape".
    print("\nStep 2: The 'peculiar shape' refers to the appearance of the coffee.")
    print("We need a coffee type whose name describes a unique mark or shape.")
    
    # Clue 3: Connecting the clues.
    print("\nStep 3: The word 'Macchiato' is Italian for 'stained' or 'marked'.")
    second_part = "Macchiato"
    print(f"The second part of the name is: {second_part}")

    # Final Answer
    print("-" * 40)
    print("A Caffè Macchiato is an espresso with a small mark, or 'stain,' of milk, which is the 'peculiar shape.'")
    print("\nThe correctly spelled name is:")
    print(f"{first_part} {second_part}")

solve_riddle()