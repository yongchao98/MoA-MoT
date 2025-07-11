import math

def calculate_iconic_ballet_turns():
    """
    This function calculates the number of turns in Kitri's most famous turning sequence.

    The role of Kitri in Don Quixote is iconic for a demanding sequence of turns.
    While the question refers to the Act I variation, the most numerically famous sequence
    is the 32 fouetté turns from the Act III grand pas de deux. This is a benchmark for ballerinas.

    We will represent the calculation for this number.
    """
    
    # The total number of turns is 32.
    total_turns = 32

    # To create an equation as requested, we can represent this number as a product.
    # For instance, we can think of it as 4 groups of 8 turns.
    groups = 4
    turns_per_group = 8

    # Calculate the total.
    calculated_total = groups * turns_per_group
    
    # Print the equation and the final answer.
    print(f"The famous sequence for Kitri consists of {calculated_total} turns.")
    print(f"The calculation can be represented by the equation: {groups} * {turns_per_group} = {calculated_total}")

calculate_iconic_ballet_turns()