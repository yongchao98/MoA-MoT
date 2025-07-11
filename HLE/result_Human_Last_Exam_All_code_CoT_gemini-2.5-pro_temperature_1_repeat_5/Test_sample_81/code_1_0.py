import sys

def analyze_tos_clauses():
    """
    Analyzes Terms of Service clauses to identify the one most likely
    to be a problematic contract of adhesion or hide material terms.
    
    This is based on a heuristic scoring system:
    - Unexpectedness (U): How surprising the term is.
    - Burdensomeness (B): How much it restricts the user.
    - Materiality (M): How much it impacts the service's core use.
    
    Problematic Score = U + B + M
    """

    clauses = {
        'A': {
            'text': "No building competitive products.",
            'u': 4, 'b': 7, 'm': 6,
            'reason': "Common in B2B SaaS, but restrictive. A user might expect this."
        },
        'B': {
            'text': "Granting a worldwide, transferable license for user content.",
            'u': 2, 'b': 8, 'm': 9,
            'reason': "Extremely broad rights grant, but now standard and expected for social media/content platforms."
        },
        'C': {
            'text': "Late charge of 'one and one-half percent (10.5%) per month'.",
            'u': 9, 'b': 8, 'm': 5,
            'reason': "Highly problematic due to the major discrepancy between the text (1.5%) and the number (10.5%), which could be deceptive. However, it only affects users who pay late."
        },
        'D': {
            'text': "No scraping or automated copying of content.",
            'u': 1, 'b': 5, 'm': 5,
            'reason': "Very standard, expected term on almost every modern website."
        },
        'E': {
            'text': "Prohibition on researching individuals in Illinois.",
            'u': 10, 'b': 9, 'm': 8,
            'reason': "Extremely specific, unusual, and material restriction that a user would not reasonably expect. It's buried in a list of more standard prohibitions and likely exists to avoid the company's legal liability under a specific state law (BIPA), directly limiting the service's function in a non-obvious way. This is a classic 'hidden material term'."
        },
        'F': {
            'text': "No using service content to train other AI models.",
            'u': 3, 'b': 7, 'm': 7,
            'reason': "A modern but increasingly standard and expected term for AI services."
        },
        'G': {
            'text': "Using user's name/profile pic in ads without compensation.",
            'u': 5, 'b': 9, 'm': 8,
            'reason': "Very burdensome and controversial, but now a well-known and expected term for major social media platforms (e.g., Facebook's 'social ads')."
        }
    }

    highest_score = -1
    best_option = None

    print("Analyzing each clause based on Unexpectedness (U), Burdensomeness (B), and Materiality (M):\n")

    for option, data in clauses.items():
        score = data['u'] + data['b'] + data['m']
        clauses[option]['score'] = score
        print(f"Option {option}: {data['text']}")
        print(f"Analysis: {data['reason']}")
        print(f"Score: {data['u']} (U) + {data['b']} (B) + {data['m']} (M) = {score}\n")

        if score > highest_score:
            highest_score = score
            best_option = option

    print("--- Conclusion ---")
    print(f"The most problematic clause is Option {best_option} with a score of {highest_score}.")
    print("It contains a highly unexpected and material term hidden within a long list of prohibitions.")
    
    # Per the instructions, print the final equation for the winning option
    winning_clause = clauses[best_option]
    print("\nFinal calculation for the winning option:")
    print(f"{winning_clause['u']} + {winning_clause['b']} + {winning_clause['m']} = {winning_clause['score']}")

if __name__ == '__main__':
    analyze_tos_clauses()