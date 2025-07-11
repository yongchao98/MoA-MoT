import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on predicted binding affinity to H2-Kd MHC.

    The ranking is determined by a scoring system based on key principles of
    MHC-peptide interaction for the H2-Kd allele:

    1.  **Anchor Residues:** Binding is dominated by anchors at Position 2 (P2) and
        Position 9 (P9).
    2.  **P2 Score:** H2-Kd prefers Tyrosine (Y) or Phenylalanine (F). Y is in the
        reference high-affinity epitope and can form a hydrogen bond that F cannot,
        so it's scored higher.
        - P2 is 'Y': +0 points (baseline)
        - P2 is 'F': -30 points (significant penalty)
    3.  **P9 Score:** H2-Kd prefers hydrophobic residues (like V) and disfavors
        charged residues (like K).
        - P9 is 'V': +0 points (baseline)
        - P9 is 'K': -50 points (very large penalty)
    4.  **Core Mutations:** Mutations from the reference high-affinity epitope (E1)
        in the core region (P3-P8) are assumed to be suboptimal. A penalty is
        applied for each mutation. A special, larger penalty is applied for Proline,
        which is known to be structurally disruptive.
        - Per mutation vs E1: -5 points
        - Additional penalty if mutation is Proline: -10 points

    The final score for each epitope starts at 100 and has penalties subtracted.
    A higher score indicates a higher predicted binding affinity.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFA LV", #Note: Space is a typo in question, will be handled
        "E5": "TYIPMRALV",
    }

    reference_seq = epitopes["E1"]
    scores = {}

    print("Calculating binding scores based on H2-Kd binding motif rules...\n")

    for name, seq in epitopes.items():
        # Sanitize sequence by removing spaces
        current_seq = seq.replace(" ", "")
        
        score = 100
        explanation = []

        # P2 anchor evaluation
        if current_seq[1] == 'F':
            score -= 30
            explanation.append("P2='F' is less optimal than 'Y' (-30)")

        # P9 anchor evaluation
        if current_seq[8] == 'K':
            score -= 50
            explanation.append("P9='K' is a highly disfavored anchor (-50)")

        # Core mutation evaluation (P3 to P8, index 2 to 7)
        if name != "E1":
            mutation_count = 0
            for i in range(2, 8):
                if current_seq[i] != reference_seq[i]:
                    mutation_count += 1
                    # Proline penalty
                    if current_seq[i] == 'P':
                        score -= 10
                        explanation.append(f"Core mutation P{i+1}='P' is structurally disruptive (-10)")
            
            mutation_penalty = mutation_count * 5
            score -= mutation_penalty
            if mutation_count > 0:
                 explanation.append(f"{mutation_count} core mutation(s) (-{mutation_penalty})")


        scores[name] = {"sequence": current_seq, "score": score, "explanation": ", ".join(explanation) or "Reference sequence"}

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1]['score'], reverse=True)

    # Print the results in a formatted table
    print("Ranking of Epitopes by Predicted Binding Affinity (High to Low):")
    print("-" * 75)
    print(f"{'Rank':<5} {'ID':<5} {'Sequence':<15} {'Score':<10} {'Reasoning':<40}")
    print("-" * 75)
    
    final_rank_order = []
    for i, (name, data) in enumerate(ranked_epitopes):
        final_rank_order.append(name)
        print(f"{i+1:<5} {name:<5} {data['sequence']:<15} {data['score']:<10} {data['explanation']:<40}")

    print("-" * 75)
    print("\nFinal Ranked Order:", ", ".join(final_rank_order))

# Execute the function
rank_epitopes()
<<<A>>>