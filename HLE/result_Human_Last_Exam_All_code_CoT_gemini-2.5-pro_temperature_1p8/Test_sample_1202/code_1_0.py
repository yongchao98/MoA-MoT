import textwrap

def analyze_gene_duplication_fates():
    """
    Analyzes the potential fates of duplicated genes to answer the multiple-choice question.
    """

    # Define the answer choices
    choices = {
        'A': 'Gene conversion',
        'B': 'Pseudogenization',
        'C': 'Neofunctionalization',
        'D': 'Subfunctionalization',
        'E': 'Adaptive radiation'
    }

    # The "equation" will be a series of logical steps leading to a conclusion.
    print("Evaluating the mechanisms for duplicate gene retention and divergence:\n")
    
    # Point 1: Analyze Gene Conversion
    point_1_num = 1
    point_1_text = f"Analysis of (A) {choices['A']}: This mechanism makes gene copies more similar. It actively works against the divergence of duplicate genes. Therefore, it is incorrect."
    
    # Point 2: Analyze Pseudogenization
    point_2_num = 2
    point_2_text = f"Analysis of (B) {choices['B']}: This is the most common fate, where one copy becomes a non-functional pseudogene. This is a mechanism of gene LOSS, not the RETENTION of two functional genes as the question asks. Therefore, it is incorrect."
    
    # Point 3: Analyze Neofunctionalization
    point_3_num = 3
    point_3_text = f"Analysis of (C) {choices['C']}: This is a key mechanism for retention. One copy keeps the original function while the other acquires a NEW beneficial function. Selection then preserves both copies, which diverge in sequence and function. This is a correct mechanism."
    
    # Point 4: Analyze Subfunctionalization
    point_4_num = 4
    point_4_text = f"Analysis of (D) {choices['D']}: This is another key mechanism. The ancestral gene's functions are partitioned between the two copies. Both are now required to fulfill the original roles, ensuring their retention and allowing them to diverge. This is also a correct mechanism."

    # Point 5: Analyze Adaptive Radiation
    point_5_num = 5
    point_5_text = f"Analysis of (E) {choices['E']}: This is a macroevolutionary pattern of species diversification, not a molecular-level mechanism for gene retention within a genome. Therefore, it is incorrect."

    # Print the "equation" of our logical steps
    print("Logical Equation leading to the answer:")
    print("-" * 35)
    print(f"Step {point_1_num}: {textwrap.fill(point_1_text, 80)}")
    print(f"Step {point_2_num}: {textwrap.fill(point_2_text, 80)}")
    print(f"Step {point_3_num}: {textwrap.fill(point_3_text, 80)}")
    print(f"Step {point_4_num}: {textwrap.fill(point_4_text, 80)}")
    print(f"Step {point_5_num}: {textwrap.fill(point_5_text, 80)}")
    print("-" * 35)
    
    # Final Conclusion from the "equation"
    final_conclusion_numbers = [point_3_num, point_4_num]
    final_conclusion_letters = ['C', 'D']
    print(f"\nCONCLUSION from steps {final_conclusion_numbers[0]} and {final_conclusion_numbers[1]}:")
    print(f"Both {choices['C']} ({final_conclusion_letters[0]}) and {choices['D']} ({final_conclusion_letters[1]}) are the primary mechanisms responsible for the retention and divergence of duplicate genes.")
    print("Both are considered correct answers in evolutionary genomics. Neofunctionalization is the classic model for the evolution of novel functions.")

analyze_gene_duplication_fates()