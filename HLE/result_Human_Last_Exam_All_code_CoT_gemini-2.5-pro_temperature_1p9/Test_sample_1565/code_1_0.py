def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """
    nodes_data = {
        'P1': {'total': 7, 'trust': 6, 'distrust': 1, 'edge_to_B': 'trust'},
        'P2': {'total': 6, 'trust': 4, 'distrust': 2, 'edge_to_B': 'trust'},
        'P3': {'total': 4, 'trust': 2, 'distrust': 2, 'edge_to_B': 'trust'},
        'N1': {'total': 6, 'trust': 3, 'distrust': 3, 'edge_to_B': 'distrust'},
        'N2': {'total': 4, 'trust': 1, 'distrust': 3, 'edge_to_B': 'distrust'},
    }

    total_score = 0
    equation_parts = []
    
    print("Calculating B's importance score step-by-step:\n")

    # Order of nodes for calculation
    node_order = ['P1', 'P2', 'P3', 'N1', 'N2']

    for node_name in node_order:
        data = nodes_data[node_name]
        total_rels = data['total']
        trust_ratings = data['trust']
        distrust_ratings = data['distrust']
        
        contribution = 0
        
        if data['edge_to_B'] == 'trust':
            # Rule 1: Positive edge
            contribution = 1 / (total_rels + 1)
            total_score += contribution
            equation_part = f"1/({total_rels}+1)"
            equation_parts.append(equation_part)
            print(f"Contribution from {node_name} (Positive Edge): {equation_part} = {contribution:.4f}")

        elif data['edge_to_B'] == 'distrust':
            # Rule 2: Negative edge with mixed ratings
            base_contribution = -1 / (total_rels + 1) * (trust_ratings / total_rels)
            
            # Rule 3: Users with more distrust than trust
            if distrust_ratings > trust_ratings:
                contribution = 1.5 * base_contribution
                total_score += contribution
                equation_part = f"1.5*(-1/({total_rels}+1) * {trust_ratings}/{total_rels})"
                equation_parts.append(equation_part)
                print(f"Contribution from {node_name} (Negative Edge, distrust > trust): {equation_part} = {contribution:.4f}")

            else:
                contribution = base_contribution
                total_score += contribution
                equation_part = f"(-1/({total_rels}+1) * {trust_ratings}/{total_rels})"
                equation_parts.append(equation_part)
                print(f"Contribution from {node_name} (Negative Edge): {equation_part} = {contribution:.4f}")

    print("\nFinal Equation:")
    final_equation_str = " + ".join(equation_parts).replace('+ -', '- ').replace('+ (', '+ (')
    print(f"B's Score = {final_equation_str}")
    
    print(f"\nFinal calculated score: {total_score:.4f}")

    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options, key=lambda k: abs(options[k] - total_score))
    print(f"The calculated score is closest to option {closest_option}: {options[closest_option]}")

if __name__ == "__main__":
    calculate_fat_score()