import operator

def evaluate_tos_clauses():
    """
    Evaluates Terms of Service clauses based on a scoring model for
    unreasonableness, reflecting concepts from contract law.
    """
    clauses = {
        'A': {'text': 'No competitive product...', 'restriction': 8, 'rights_grab': 2, 'surprise': 4},
        'B': {'text': 'You grant us a worldwide license to your content...', 'restriction': 3, 'rights_grab': 7, 'surprise': 5},
        'C': {'text': 'Late charges and non-refundable fees...', 'restriction': 2, 'rights_grab': 1, 'surprise': 1},
        'D': {'text': 'No scraping or monitoring...', 'restriction': 5, 'rights_grab': 1, 'surprise': 2},
        'E': {'text': 'List of prohibitions including no Illinois research...', 'restriction': 7, 'rights_grab': 1, 'surprise': 9},
        'F': {'text': 'No abuse, disruption, or AI training...', 'restriction': 4, 'rights_grab': 1, 'surprise': 3},
        'G': {'text': 'Use of your identity in ads without compensation...', 'restriction': 5, 'rights_grab': 9, 'surprise': 8}
    }

    # Weights for the scoring formula
    w_restriction = 1.0
    w_rights = 2.0
    w_surprise = 1.5

    results = {}

    print("Evaluating TOS clauses based on a quantitative model...\n")

    for key, values in clauses.items():
        restriction = values['restriction']
        rights_grab = values['rights_grab']
        surprise = values['surprise']

        # Calculate the total unreasonableness score
        total_score = (restriction * w_restriction) + (rights_grab * w_rights) + (surprise * w_surprise)
        results[key] = total_score

        # Print the equation for each calculation
        print(f"Clause {key}:")
        print(f"  Calculation: ({restriction} * {w_restriction}) + ({rights_grab} * {w_rights}) + ({surprise} * {w_surprise}) = {total_score}")
        print("-" * 20)

    # Find the clause with the highest score
    most_likely_candidate = max(results.items(), key=operator.itemgetter(1))

    print(f"\nConclusion:")
    print(f"The clause with the highest unreasonableness score is '{most_likely_candidate[0]}' with a score of {most_likely_candidate[1]}.")
    print("This clause permits the use of a user's identity (username, profile picture, actions) in connection with advertising without any compensation, which is a significant appropriation of rights that a user may not reasonably expect.")

evaluate_tos_clauses()
<<<G>>>