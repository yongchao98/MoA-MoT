import textwrap

def analyze_immune_receptors():
    """
    Analyzes potential receptors to engineer a T cell into an
    Antigen-Presenting Cell (APC) and prints the reasoning.
    """
    # Define the biological problem
    problem_statement = (
        "The goal is to engineer a T cell to act as an Antigen-Presenting Cell (APC). "
        "A key function of an APC is to process external antigens and present them to other T cells. "
        "We need to identify which receptor would enable this function."
    )

    # Dictionary of options with their biological functions
    options = {
        'A': {
            'name': 'CD86',
            'function': 'A co-stimulatory molecule found on APCs. It provides a crucial second signal to T cells for full activation, but it does not present the antigen.'
        },
        'B': {
            'name': 'CD80',
            'function': 'Another co-stimulatory molecule on APCs, functionally similar to CD86. It is necessary for robust T cell activation but is not the antigen-presenting molecule.'
        },
        'C': {
            'name': 'MHC class I',
            'function': 'Presents internal (endogenous) peptides, such as viral proteins, to cytotoxic T cells. T cells already express MHC class I, so adding more wouldn\'t grant them the ability to present external antigens like a professional APC.'
        },
        'D': {
            'name': 'TIM-4',
            'function': 'A receptor primarily involved in the phagocytosis (engulfment) of apoptotic (dying) cells. It is not directly involved in presenting antigens to T cells.'
        },
        'E': {
            'name': 'MHC class II',
            'function': 'The quintessential molecule of professional APCs (like macrophages and dendritic cells). It specifically binds to peptides from external (exogenous) sources and presents them to helper T cells, initiating a broader adaptive immune response.'
        }
    }

    # Print the step-by-step analysis
    print("Step-by-Step Analysis:")
    print("="*60)
    for line in textwrap.wrap(problem_statement, width=70):
        print(line)
    print("\nEvaluating the options:")
    print("-"*60)

    for key, value in options.items():
        print(f"Option {key} ({value['name']}):")
        for line in textwrap.wrap(value['function'], width=65, initial_indent='  ', subsequent_indent='  '):
            print(line)
        print("")

    # Print the conclusion
    conclusion = (
        "To enable a T cell to act as an APC, it must be given the molecular machinery "
        "to present antigens it has taken up from the environment. MHC class II is the specific "
        "receptor responsible for this function. Therefore, engineering a T cell to express "
        "MHC class II would allow it to present external antigens to other T cells."
    )
    print("Conclusion:")
    print("="*60)
    for line in textwrap.wrap(conclusion, width=70):
        print(line)

    print("\nFinal Answer Choice: E")

if __name__ == '__main__':
    analyze_immune_receptors()