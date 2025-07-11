def analyze_bansenshukai_pattern():
    """
    This function analyzes the symbolic pattern from the Bansenshukai scroll
    and presents the numerical representation as an equation.
    """
    # The pattern is ⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤.
    # We can interpret this as sequential groups of same-colored circles.
    # Group 1: ⬤ (1 black circle)
    # Group 2: ○○ (2 white circles)
    # Group 3: ⬤⬤⬤⬤ (4 black circles)
    # Group 4: ○ (1 white circle)
    # Group 5: ⬤⬤⬤⬤⬤ (5 black circles)
    
    # The numbers derived from the pattern's groups.
    numbers = [1, 2, 4, 1, 5]
    
    # The total number of symbols in the pattern.
    total_symbols = sum(numbers)
    
    # The prompt requires outputting each number in the final equation.
    # We will format this as a sum.
    equation_str = " + ".join(map(str, numbers))
    
    print(f"The symbolic pattern contains {total_symbols} circles in total.")
    print("The numerical sequence derived from the groups of circles is represented by the equation:")
    print(f"{equation_str} = {total_symbols}")

analyze_bansenshukai_pattern()