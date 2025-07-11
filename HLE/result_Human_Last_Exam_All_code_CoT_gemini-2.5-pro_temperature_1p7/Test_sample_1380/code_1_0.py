def solve_bachelier_question():
    """
    This script explains the physical process Louis Bachelier's work described.
    """
    
    # Year of Louis Bachelier's thesis "Théorie de la Spéculation"
    year_of_thesis = 1900
    
    # The name of the physical process
    process_name = "Brownian motion"
    
    # Explanation constructing the final answer
    explanation = (
        f"In his thesis published in the year {year_of_thesis}, Louis Bachelier demonstrated that the probabilities\n"
        f"of stock price fluctuations spread over time in a manner that is mathematically identical\n"
        f"to a specific physical process.\n\n"
        f"This process is the random, erratic movement of microscopic particles suspended in a fluid,\n"
        f"and it is named: {process_name}.\n\n"
        f"Bachelier's work was the first rigorous mathematical model of this phenomenon, predating\n"
        f"Albert Einstein's famous paper on the subject by five years."
    )
    
    print(explanation)

solve_bachelier_question()