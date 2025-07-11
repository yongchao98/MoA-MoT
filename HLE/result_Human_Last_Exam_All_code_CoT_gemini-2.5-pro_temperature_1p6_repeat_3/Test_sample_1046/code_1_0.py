import sys

def analyze_bacterial_resistance():
    """
    Analyzes different scenarios for acquiring drug resistance in bacteria
    to determine which one can match the pace of lateral gene transfer (LGT).
    """

    print("Problem: A bacterium with a stable genome (Bacterium 2) acquires drug resistance as fast as a bacterium with frequent lateral transfer (Bacterium 1). How?")
    print("-" * 80)
    print("Let's analyze the options by scoring the evolutionary mechanisms for their speed and impact.")
    print("\nThe goal is to find a mechanism for Bacterium 2 that is powerful enough to match the speed of LGT.\n")

    # Define scores for different evolutionary events.
    # LGT is the benchmark for speed. We need a combination of factors to match it.
    scores = {
        'initial_mutation': 1,      # Essential, but slow on its own
        'compensation': 3,          # Recovers fitness, crucial for spread
        'high_fitness_gain': 5,     # The main accelerator, allows rapid population takeover
        'cross_resistance': 4,      # Very advantageous, accelerates selection
        'contamination_theory': -10 # Dismisses the biological question
    }

    # Define the answer choices and their components
    options = {
        'A': {
            'description': "Rare mutations occurred causing resistance.",
            'components': ['initial_mutation']
        },
        'B': {
            'description': "Acquired compensatory mutations that greatly increased fitness and led to cross-resistance.",
            'components': ['initial_mutation', 'compensation', 'high_fitness_gain', 'cross_resistance']
        },
        'C': {
            'description': "Contamination, since LGT is the major way to acquire resistance.",
            'components': ['contamination_theory']
        },
        'D': {
            'description': "Mutations led to cross-resistance but without compensatory mutations.",
            'components': ['initial_mutation', 'cross_resistance']
        },
        'E': {
            'description': "Acquired compensatory mutations that followed rare resistance mutations.",
            'components': ['initial_mutation', 'compensation']
        }
    }

    best_option = ''
    max_score = -sys.maxsize

    # Analyze each option
    for option_key, data in options.items():
        print(f"--- Analyzing Option {option_key} ---")
        print(f"Description: {data['description']}")
        
        current_score = 0
        equation_parts = []
        for component in data['components']:
            score = scores.get(component, 0)
            current_score += score
            equation_parts.append(str(score))
        
        # Displaying the final equation with each number
        equation_str = " + ".join(equation_parts)
        print(f"Speed Score Calculation: {equation_str} = {current_score}")
        
        if option_key == 'A':
            print("Analysis: A rare mutation is the starting point, but it's typically slow to spread due to fitness costs. This doesn't explain the rapid pace.\n")
        elif option_key == 'B':
            print("Analysis: This option combines multiple powerful factors. A high fitness gain allows the new strain to rapidly outcompete others, and cross-resistance provides a major selective advantage. This combination is a strong candidate for matching the speed of LGT.\n")
        elif option_key == 'C':
            print("Analysis: This is an external explanation, not a biological one. It avoids answering the core evolutionary question.\n")
        elif option_key == 'D':
            print("Analysis: Cross-resistance is beneficial, but without compensation for the likely fitness cost of the initial mutation, the spread would be hindered, not accelerated.\n")
        elif option_key == 'E':
            print("Analysis: Compensation is important for preventing the resistant strain from dying out, but without a major fitness gain, it wouldn't explain an exceptionally *fast* spread.\n")
            
        if current_score > max_score:
            max_score = current_score
            best_option = option_key

    print("-" * 80)
    print(f"Conclusion: Option {best_option} has the highest score ({max_score}).")
    print("It provides the most comprehensive mechanism (initial mutation, fitness compensation, significant fitness enhancement, and cross-resistance) to explain how a population with a stable genome can achieve rapid resistance spread, rivaling the pace of lateral gene transfer.")

analyze_bacterial_resistance()