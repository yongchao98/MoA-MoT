def analyze_genetic_alterations():
    """
    Analyzes genetic alterations to identify the one matching
    the criteria for Copy-Neutral Loss of Heterozygosity (CN-LOH).
    """
    alterations = {
        'A': {
            'name': 'Mitotic recombination',
            'copy_number_change': 'neutral',
            'causes_LOH': True,
            'note': 'A major mechanism causing CN-LOH.'
        },
        'B': {
            'name': 'A deletion of a chromosomal region',
            'copy_number_change': 'loss',
            'causes_LOH': True,
            'note': 'Not copy-neutral.'
        },
        'C': {
            'name': 'Trisomy',
            'copy_number_change': 'gain',
            'causes_LOH': False,
            'note': 'Not copy-neutral.'
        },
        'D': {
            'name': 'Uniparental disomy',
            'copy_number_change': 'neutral',
            'causes_LOH': True,
            'note': 'The definitional state of CN-LOH.'
        },
        'E': {
            'name': 'Duplication of a chromosomal region',
            'copy_number_change': 'gain',
            'causes_LOH': False,
            'note': 'Not copy-neutral.'
        }
    }

    print("Analyzing options for Copy-Neutral Loss of Heterozygosity (CN-LOH):")
    print("CN-LOH requires two conditions: 1. Copy number is neutral (remains 2). 2. Heterozygosity is lost.")
    print("-" * 75)

    best_choice = None
    best_explanation = ""

    for choice, props in alterations.items():
        is_cn_loh = props['copy_number_change'] == 'neutral' and props['causes_LOH']
        print(f"Choice {choice}: {props['name']}")
        print(f"  - Is it copy-neutral? {'Yes' if props['copy_number_change'] == 'neutral' else 'No'} (It's a copy number {props['copy_number_change']})")
        print(f"  - Does it cause LOH? {'Yes' if props['causes_LOH'] else 'No'}")
        print(f"  - Meets criteria for CN-LOH? {'Yes' if is_cn_loh else 'No'}")
        print(f"  - Note: {props['note']}")
        print("-" * 75)
        
        # Uniparental Disomy is the most fundamental definition of CN-LOH
        if props['name'] == 'Uniparental disomy':
            best_choice = choice
            best_explanation = props['note']
            
    print("\nConclusion:")
    print("While Mitotic Recombination is a mechanism for CN-LOH, Uniparental Disomy (UPD)")
    print("is the term that describes the resulting state itself: having two copies (disomy)")
    print("of a chromosome from a single parent (uniparental). It is the most direct and")
    print("accurate answer describing copy-neutral loss of heterozygosity.")
    print(f"\nThe best choice is '{best_choice}' because it represents {best_explanation}")

if __name__ == '__main__':
    analyze_genetic_alterations()