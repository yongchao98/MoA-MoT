import sys

def solve_paleoclimate_question():
    """
    Analyzes climate factors to determine which combination leads to low 18O isotope levels.
    """
    # Define the answer choices with their corresponding factors
    options = {
        'A': ('Wet', 'Warm', 'Shallow lake levels'),
        'B': ('Dry', 'Warm', 'Shallow lake levels'),
        'C': ('Wet', 'Cold', 'Shallow lake levels'),
        'D': ('Dry', 'Cold', 'Shallow lake levels'),
        'E': ('Wet', 'Warm', 'High lake levels'),
        'F': ('Dry', 'Warm', 'High lake levels'),
        'G': ('Wet', 'Cold', 'High lake levels'),
        'H': ('Dry', 'Cold', 'High lake levels')
    }

    # Define scores for each factor's effect on 18O levels.
    # Positive score = promotes LOW 18O
    # Negative score = promotes HIGH 18O
    factor_scores = {
        'Wet': 1,
        'Dry': -1,
        'Cold': 1,
        'Warm': -1,
        'High lake levels': 1,
        'Shallow lake levels': -1
    }

    print("Analyzing factors for low 18O levels in sulfate deposits.")
    print("Scoring: +1 for factors promoting low 18O, -1 for factors promoting high 18O.\n")

    best_option = None
    max_score = -sys.maxsize

    # Calculate the score for each option
    for option, factors in options.items():
        score = sum(factor_scores[f] for f in factors)
        if score > max_score:
            max_score = score
            best_option = option

    # Print the analysis of the best option
    winning_factors = options[best_option]
    factor1, factor2, factor3 = winning_factors
    score1 = factor_scores[factor1]
    score2 = factor_scores[factor2]
    score3 = factor_scores[factor3]
    
    print(f"The best combination of factors for low 18O is Option {best_option}:")
    print(f"'{factor1}, {factor2}, {factor3}'")
    
    # Display the "equation" as requested
    print("\nScore Calculation:")
    print(f"{factor1} ({score1}) + {factor2} ({score2}) + {factor3} ({score3}) = {max_score}")
    
    print("\nThis combination minimizes evaporation (Cold, Wet) and maximizes inflow of 18O-depleted water (Wet, High lake levels), leading to low 18O in the lake.")

solve_paleoclimate_question()