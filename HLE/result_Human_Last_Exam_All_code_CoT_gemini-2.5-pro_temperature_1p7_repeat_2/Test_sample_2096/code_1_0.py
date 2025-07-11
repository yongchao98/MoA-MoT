def solve_riddle_and_calculate():
    """
    This function explains the riddle's logic and performs a calculation
    based on the years mentioned in the problem description.
    """
    # The key years from the riddle
    year_written_clue = 1960
    year_constantinople = 1453

    # The riddle's answer is a word, not a number.
    # The word is "ILLITERATE".
    # Pope Paul II was accused by his enemies of being hostile to learning.
    # To a Renaissance leader, being illiterate would be a source of great shame.
    # The clue "written in the 1960s" refers to the fact that in the USSR
    # during that decade, illiteracy was a major social issue that was heavily
    # written about in propaganda and public campaigns, making it a "hot topic".

    # Per the instructions, create and print an equation using the numbers.
    difference = year_written_clue - year_constantinople

    print("The riddle combines historical facts with cultural context.")
    print("Here is a simple calculation using the two key dates from the text:")
    
    # Printing each part of the equation as requested
    print(year_written_clue, end=' ')
    print('-', end=' ')
    print(year_constantinople, end=' ')
    print('=', end=' ')
    print(difference)

solve_riddle_and_calculate()