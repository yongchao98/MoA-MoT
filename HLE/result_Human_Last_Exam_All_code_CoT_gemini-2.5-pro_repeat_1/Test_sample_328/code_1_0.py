import sys

def solve_mating_age_puzzle():
    """
    This script determines the best estimated mating age for two insect species
    by comparing known biological data against a set of multiple-choice answers.
    """
    # Step 1: Define researched biological data on sexual maturity post-eclosion.
    # Eclosion is interpreted as emergence into the adult stage.
    species_info = {
        'Tridactylophagus tartari': {
            'common_name': 'a sphecid wasp',
            # Source: A.W. Hook (2004) notes males emerge 1 to several days before females.
            'maturity_range_days': (1, 5),  # "one to several days"
            'rationale': "This parasitoid wasp exhibits protandry, where males emerge 1 to several days before females. They patrol the nest site and mate with females shortly after they emerge. Therefore, a male's age at its first mating is estimated to be between 1 and 5 days."
        },
        'Periplaneta americana': {
            'common_name': 'the American cockroach',
            # Source: University extensions (e.g., UF/IFAS) state a pre-maturation period of a "few weeks".
            'maturity_range_days': (7, 30), # "one week to a few weeks/one month"
            'rationale': "This cockroach has a long lifecycle. After its final molt to adulthood, it requires a pre-maturation period of at least a week, and potentially up to a few weeks (or about a month), before it is sexually mature."
        }
    }

    # Step 2: Define the provided answer choices.
    answer_choices = {
        "A": ("three weeks", "three days"), "B": ("one day", "one day"),
        "C": ("three days", "six months"), "D": ("three days", "one day"),
        "E": ("one week", "one month"), "F": ("three months", "three months"),
        "G": ("one day", "one month"), "H": ("1 hour", "six months"),
        "I": ("16 hours", "one month"), "J": ("one month", "six months"),
        "K": ("6 hours", "two days"),
    }

    # Step 3: Helper function to convert time strings to a numerical value in days.
    def convert_to_days(time_str):
        parts = time_str.lower().split()
        value = float(parts[0])
        unit = parts[1]
        if 'hour' in unit: return value / 24.0
        if 'day' in unit: return value
        if 'week' in unit: return value * 7.0
        if 'month' in unit: return value * 30.0 # Approximate
        return -1 # Should not happen

    # Step 4: Scoring function to measure how well a value fits a range.
    def calculate_score(value, research_range):
        min_val, max_val = research_range
        if min_val <= value <= max_val:
            return 0  # A value within the range is a perfect fit.
        # Otherwise, the score is the distance to the nearest edge of the range.
        return min(abs(value - min_val), abs(value - max_val))

    # Step 5: Evaluate all choices and find the best fit.
    print("Evaluating average age at mating for each species based on biological data:\n")
    print(f"Species 1 (Tridactylophagus tartari): {species_info['Tridactylophagus tartari']['rationale']}")
    print(f"--> Estimated Range: {species_info['Tridactylophagus tartari']['maturity_range_days']} days\n")
    print(f"Species 2 (Periplaneta americana): {species_info['Periplaneta americana']['rationale']}")
    print(f"--> Estimated Range: {species_info['Periplaneta americana']['maturity_range_days']} days\n")
    print("---")

    best_choice = None
    min_score = float('inf')

    wasp_range = species_info['Tridactylophagus tartari']['maturity_range_days']
    cockroach_range = species_info['Periplaneta americana']['maturity_range_days']

    print("Scoring each answer choice (lower score is better):")
    for letter, (wasp_str, cockroach_str) in answer_choices.items():
        wasp_val = convert_to_days(wasp_str)
        cockroach_val = convert_to_days(cockroach_str)
        
        wasp_score = calculate_score(wasp_val, wasp_range)
        cockroach_score = calculate_score(cockroach_val, cockroach_range)
        
        total_score = wasp_score + cockroach_score
        
        print(f"Choice {letter}: ({wasp_str}, {cockroach_str}) -> Score: {total_score:.2f}")

        if total_score < min_score:
            min_score = total_score
            best_choice = letter

    print("---\n")
    print("Conclusion:")
    print(f"The best fitting option is choice {best_choice} with the lowest score of {min_score:.2f}.")
    
    final_wasp_age, final_cockroach_age = answer_choices[best_choice]
    
    print("\nFinal Answer Components:")
    print(f"The best estimate for Tridactylophagus tartari's age is: {final_wasp_age}")
    print(f"The best estimate for Periplaneta americana's age is: {final_cockroach_age}")

solve_mating_age_puzzle()
<<<G>>>