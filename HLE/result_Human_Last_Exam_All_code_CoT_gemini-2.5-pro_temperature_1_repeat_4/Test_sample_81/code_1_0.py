import operator

def analyze_tos_clauses():
    """
    Analyzes terms of service clauses to identify the one most likely to be
    part of a contract of adhesion or hide a material term by being highly
    unexpected.
    """
    
    clauses = {
        'A': {
            'description': "No-compete / No-benchmark clause.",
            'surprise_factor': 3,
            'reasoning': "Common in B2B agreements, but can be restrictive. Somewhat expected."
        },
        'B': {
            'description': "Broad license grant for user-submitted content.",
            'surprise_factor': 2,
            'reasoning': "Standard and necessary for social media/content platforms to function. Widely expected."
        },
        'C': {
            'description': "Late fees and non-refundable payments.",
            'surprise_factor': 1,
            'reasoning': "Very standard financial term in service agreements. Highly expected."
        },
        'D': {
            'description': "Prohibition on scraping or automated data copying.",
            'surprise_factor': 1,
            'reasoning': "Extremely common and standard clause to protect a service's data. Highly expected."
        },
        'E': {
            'description': "List of prohibitions, including a ban on researching Illinois residents.",
            'surprise_factor': 10,
            'reasoning': ("The prohibition against researching individuals in a single specific state (Illinois) is "
                          "extremely unusual and unexpected. A user would have no reason to anticipate this. "
                          "It points to the company avoiding a specific, significant legal liability (likely BIPA) "
                          "in a way that fails the 'reasonable expectations' test for contracts of adhesion.")
        },
        'F': {
            'description': "Prohibition on using content to train AI models.",
            'surprise_factor': 4,
            'reasoning': "A modern, but increasingly common and expected clause for AI and data-rich services."
        },
        'G': {
            'description': "Permission to use user profile/actions in connection with ads.",
            'surprise_factor': 6,
            'reasoning': ("The business model for many 'free' social networks. While a major term, it is "
                          "unfortunately becoming a reasonably expected one for such services.")
        }
    }

    print("Analyzing contract clauses based on the 'doctrine of reasonable expectations'.")
    print("A higher 'surprise factor' means the term is less likely to be reasonably expected by a user.\n")

    for key, value in clauses.items():
        print(f"Clause {key}: {value['description']}")
        print(f"Surprise Factor: {value['surprise_factor']}/10")
        print(f"Reasoning: {value['reasoning']}\n")

    # Sort clauses by surprise factor in descending order
    sorted_clauses = sorted(clauses.items(), key=lambda item: item[1]['surprise_factor'], reverse=True)
    
    print("--- Conclusion ---")
    print("Clause E has the highest surprise factor because a geographic restriction on researching people in Illinois is not a term a reasonable person would ever expect in a standard contract.")
    
    # Create and print the "equation" as requested
    equation_parts = [str(item[1]['surprise_factor']) for item in sorted_clauses]
    print("\nFinal Equation (Surprise Factors Ranked):")
    # This prints each number in the final equation
    for i, part in enumerate(equation_parts):
        print(part, end="")
        if i < len(equation_parts) - 1:
            print(" > ", end="")
    print("\n")


analyze_tos_clauses()

# The final answer is the letter corresponding to the highest surprise factor.
print("<<<E>>>")