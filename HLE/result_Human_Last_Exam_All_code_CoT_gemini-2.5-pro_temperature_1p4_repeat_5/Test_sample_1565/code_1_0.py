def calculate_fat_score():
    """
    Calculates and prints the FAT score for user B based on the provided data and rules.
    """
    # Data for users connecting to B
    # Format: {name: [total_rels, trust_ratings, distrust_ratings, edge_type_to_B]}
    users = {
        'P1': [7, 6, 1, 'positive'],
        'P2': [6, 4, 2, 'positive'],
        'P3': [4, 2, 2, 'positive'],
        'N1': [6, 3, 3, 'negative'],
        'N2': [4, 1, 3, 'negative']
    }

    total_score = 0
    contributions = {}

    print("Calculating B's FAT Score:")
    print("--------------------------------")

    # Calculate contribution for each user
    for name, data in users.items():
        total_rels, trust_ratings, distrust_ratings, edge_type = data
        contribution = 0

        if edge_type == 'positive':
            # Rule 1: Positive edge contributes 1/(total_relationships + 1)
            contribution = 1 / (total_rels + 1)
            print(f"Contribution from {name} (Positive Edge):")
            print(f"  = 1 / ({total_rels} + 1)")
            print(f"  = {contribution:.4f}\n")
        
        elif edge_type == 'negative':
            # Rule 2: Negative edge with mixed ratings
            contribution = -1 / (total_rels + 1) * (trust_ratings / total_rels)
            
            print(f"Contribution from {name} (Negative Edge):")
            print(f"  Base = -1 / ({total_rels} + 1) * ({trust_ratings} / {total_rels})")

            # Rule 3: Users with more distrust than trust get 1.5x negative weight
            if distrust_ratings > trust_ratings:
                print(f"  User {name} has more distrust ({distrust_ratings}) than trust ({trust_ratings}), applying 1.5x weight.")
                contribution *= 1.5
                print(f"  Final = ({contribution/1.5:.4f}) * 1.5")

            print(f"  = {contribution:.4f}\n")

        contributions[name] = contribution
        total_score += contribution

    # Print the final equation
    print("--------------------------------")
    print("Final Score Calculation:")
    equation_parts = []
    for name, value in contributions.items():
        if value < 0:
            equation_parts.append(f"({value:.4f})")
        else:
            equation_parts.append(f"{value:.4f}")
    
    equation_str = " + ".join(equation_parts)
    print(f"Score = {equation_str}")
    print(f"Total FAT Score for B = {total_score:.4f}")

# Run the calculation
calculate_fat_score()