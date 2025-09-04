def check_chemistry_answer(selected_option: str):
    """
    This function checks the correctness of the selected option for a multi-step
    organic chemistry synthesis problem.

    The question is to identify the number of 13C-NMR signals for the final product E.
    """
    # The final product E is 3,4-diethylhex-3-ene, with the structure (CH3CH2)2C=C(CH2CH3)2.
    # Let's analyze the symmetry to determine the number of unique carbon environments.
    #
    # Structure:
    #       Et
    #      /
    #   Et-C=C-Et
    #      \
    #       Et
    #
    # The molecule is highly symmetrical.
    # 1. The two carbons of the central double bond (>C=C<) are equivalent. (Signal 1)
    # 2. All four ethyl groups are chemically equivalent due to the symmetry.
    # 3. Within each ethyl group, the methylene (-CH2-) and methyl (-CH3) carbons are distinct.
    # 4. Therefore, all four methylene carbons are equivalent to each other. (Signal 2)
    # 5. And all four methyl carbons are equivalent to each other. (Signal 3)
    #
    # Total number of signals = 3.
    correct_number_of_signals = 3

    # The options provided in the question.
    options = {
        'A': 3,
        'B': 11,
        'C': 6,
        'D': 8
    }

    # Check if the selected option is valid.
    if selected_option.upper() not in options:
        return f"Invalid option '{selected_option}'. Please choose from A, B, C, or D."

    # Retrieve the value for the selected option.
    submitted_value = options[selected_option.upper()]

    # Compare the submitted value with the correct answer.
    if submitted_value == correct_number_of_signals:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The selected option corresponds to {submitted_value} signals.\n"
            f"The final product of the reaction sequence is 3,4-diethylhex-3-ene.\n"
            f"Due to the high symmetry of this molecule, there are only 3 unique carbon environments:\n"
            f"1. The two equivalent carbons of the C=C double bond.\n"
            f"2. The four equivalent -CH2- carbons of the ethyl groups.\n"
            f"3. The four equivalent -CH3 carbons of the ethyl groups.\n"
            f"Therefore, the correct number of 13C-NMR signals is 3."
        )
        return reason

# Example of how to use the checker:
# To check if option 'A' is correct:
# print(check_chemistry_answer('A'))
# >>> Correct

# To check if option 'C' is correct:
# print(check_chemistry_answer('C'))
# >>> Incorrect. The selected option corresponds to 6 signals.
#     The final product of the reaction sequence is 3,4-diethylhex-3-ene.
#     Due to the high symmetry of this molecule, there are only 3 unique carbon environments:
#     1. The two equivalent carbons of the C=C double bond.
#     2. The four equivalent -CH2- carbons of the ethyl groups.
#     3. The four equivalent -CH3 carbons of the ethyl groups.
#     Therefore, the correct number of 13C-NMR signals is 3.