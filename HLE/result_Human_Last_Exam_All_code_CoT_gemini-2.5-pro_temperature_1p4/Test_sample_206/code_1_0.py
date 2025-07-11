def format_creature_list(creature_numbers):
    """Formats a list of numbers into the required string format '(n1), (n2), ...'."""
    # Ensure numbers are sorted in increasing order
    creature_numbers.sort()
    
    # Create the formatted string for each number
    formatted_parts = [f"({num})" for num in creature_numbers]
    
    # Join the parts with a comma and space
    return ", ".join(formatted_parts)

def solve_mtg_puzzle():
    """
    Solves the MTG puzzle by identifying the optimal attack and the resulting casualties.
    """
    # Based on the analysis, the optimal play is to attack with all creatures.
    attacking_creatures = [1, 2, 3, 4]

    # When played optimally, the opponent blocks to survive but loses all their creatures.
    # Two of our creatures also die in the exchange.
    # My creatures that die: (2) Centaur Courser, (4) Wind Drake
    # Opponent's creatures that die: (5) River Bear, (6) Grizzly Bears, (7) Rusted Sentinel, (8) Skywinder Drake
    dead_creatures = [2, 4, 5, 6, 7, 8]

    # Format both lists into strings
    attackers_str = format_creature_list(attacking_creatures)
    dead_creatures_str = format_creature_list(dead_creatures)

    # Print the final answer in the specified format
    print(f"{attackers_str}; {dead_creatures_str}")

solve_mtg_puzzle()