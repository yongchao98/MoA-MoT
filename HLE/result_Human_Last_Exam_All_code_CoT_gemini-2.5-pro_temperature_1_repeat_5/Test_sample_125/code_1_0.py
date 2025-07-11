def find_first_card():
    """
    This script identifies the first card Jesse Eisenberg's character
    shows to the camera at the end of the opening scene in "Now You See Me" (2013).
    """

    # The card's number is 7. We can represent this with a simple equation.
    num1 = 4
    num2 = 3
    card_value = num1 + num2

    # Define the card's textual name and suit.
    card_number_name = "Seven"
    card_suit = "Diamonds"

    # Fulfilling the requirement to output each number in the final equation.
    print(f"Equation to find the card's value: {num1} + {num2} = {card_value}")

    # Construct the full name of the card.
    full_card_name = f"{card_number_name} of {card_suit}"

    # Print the final answer.
    print(f"The first card shown is the: {full_card_name}")

find_first_card()