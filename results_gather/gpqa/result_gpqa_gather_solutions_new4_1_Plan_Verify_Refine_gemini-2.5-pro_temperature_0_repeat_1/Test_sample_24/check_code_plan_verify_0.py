def check_correctness():
    # Define the options from the question
    options = {
        'A': {
            'A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'B': '4-methyl-1-phenylpent-3-en-1-one'
        },
        'B': {
            'A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        },
        'C': {
            'A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'B': '4-methyl-1-phenylpent-3-en-1-one'
        },
        'D': {
            'A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        }
    }

    # The answer to check
    llm_answer = 'B'

    # --- Step 1: Analyze Reaction A ---
    # Reaction: A + H2SO4 ---> 2,8-dimethylspiro[4.5]decan-6-one
    # This is a Pinacol rearrangement, which requires a 1,2-diol as the starting material.
    # The name 'diol' indicates two alcohol groups. The alternative is a single 'ol'.
    correct_reactant_A_name = '2,7-dimethyloctahydronaphthalene-4a,8a-diol'
    
    # Check if the reactant A in the selected option is correct.
    selected_reactant_A = options[llm_answer]['A']
    if selected_reactant_A != correct_reactant_A_name:
        reason = (
            f"Incorrect reactant for A. The selected option '{llm_answer}' proposes "
            f"A = '{selected_reactant_A}'.\n"
            "The reaction A + H2SO4 ---> 2,8-dimethylspiro[4.5]decan-6-one is a Pinacol rearrangement, "
            "which requires a 1,2-diol as the starting material. "
            f"The correct reactant is '{correct_reactant_A_name}'. "
            f"The proposed reactant '{selected_reactant_A}' is a secondary alcohol, which would undergo dehydration, not rearrangement to a ketone."
        )
        return reason

    # --- Step 2: Analyze Reaction B ---
    # Reaction: B + BuLi + H+ ---> 4-methyl-1-phenylpent-3-en-1-ol
    # This is a Wittig rearrangement, which starts from an ether. BuLi acts as a base.
    # The product does not contain a butyl group, ruling out nucleophilic addition to a ketone.
    # The name 'oxy' indicates an ether linkage. The alternative is a 'one' (ketone).
    correct_reactant_B_name = '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'

    # Check if the reactant B in the selected option is correct.
    selected_reactant_B = options[llm_answer]['B']
    if selected_reactant_B != correct_reactant_B_name:
        reason = (
            f"Incorrect reactant for B. The selected option '{llm_answer}' proposes "
            f"B = '{selected_reactant_B}'.\n"
            "The reaction B + BuLi + H+ ---> 4-methyl-1-phenylpent-3-en-1-ol is a Wittig rearrangement, "
            "which requires an ether as the starting material. "
            f"The correct reactant is '{correct_reactant_B_name}'. "
            f"The proposed reactant '{selected_reactant_B}' is a ketone. Reacting it with BuLi would cause nucleophilic addition of a butyl group, "
            "leading to a different product (a tertiary alcohol with a butyl group)."
        )
        return reason

    # --- Step 3: Final Conclusion ---
    # If both reactants in the selected option are correct, the answer is correct.
    return "Correct"

# Execute the check
# print(check_correctness())