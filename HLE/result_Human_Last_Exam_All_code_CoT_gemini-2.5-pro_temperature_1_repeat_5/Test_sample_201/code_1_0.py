def solve_mosquito_problem():
    """
    Calculates a relative abundance score for different ponds to identify the one
    posing the greatest medical threat.
    """
    # Ponds are defined by their characteristics from the answer choices.
    # Format: {Choice: (size_in_feet, age_in_years)}
    ponds = {
        'A': (10, 1),
        'C': (30, 1),
        'D': (10, 5),
        'E': (30, 5)
    }

    print("To find the pond with the highest mosquito abundance, we'll score each pond.")
    print("The model assumes abundance increases with both size and age.")
    print("Abundance Score = Pond Size (ft) * Pond Age (years)\n")

    best_pond = None
    max_score = -1
    
    # We iterate through the choices in order to present them clearly.
    for choice in sorted(ponds.keys()):
        size, age = ponds[choice]
        # The core calculation for the abundance score
        score = size * age
        
        print(f"Pond {choice} ({size} feet square, {age} year(s) old):")
        # The final equation with numbers, as requested
        print(f"Score = {size} * {age} = {score}")

        if score > max_score:
            max_score = score
            best_pond = choice

    print(f"\nConclusion: Pond {best_pond} has the highest abundance score.")
    print("Therefore, the 30 feet square, five-year-old pond represents the greatest medical threat.")

solve_mosquito_problem()