def check_pinacol_rearrangement_answer():
    """
    Checks the correctness of the selected answer for the Pinacol rearrangement question.
    The function models the reaction rules for carbocation stability and migratory aptitude.
    """

    # The final answer provided is 'C'. Let's define the products from option C.
    answer_to_check = {
        'A': '3-ethyl-3-phenylpentan-2-one',
        'B': '3-(4-hydroxyphenyl)-3-phenylpentan-2-one',
        'C': '2,2,2-tris(4-methoxyphenyl)-1-phenylethan-1-one'
    }

    # --- Define chemical groups with properties for stability and migration ---
    # 'stab_score': Contribution to carbocation stability.
    # 'mig_apt': Migratory aptitude rank.
    groups = {
        'Me':    {'name': 'methyl',           'stab_score': 0, 'mig_apt': 1},
        'Et':    {'name': 'ethyl',            'stab_score': 0, 'mig_apt': 2},
        'Ph':    {'name': 'phenyl',           'stab_score': 10, 'mig_apt': 10}, # Benzylic
        'ArOH':  {'name': '4-hydroxyphenyl',  'stab_score': 15, 'mig_apt': 11}, # Benzylic + EDG
        'ArOMe': {'name': '4-methoxyphenyl',  'stab_score': 15, 'mig_apt': 12}, # Benzylic + EDG
    }

    # --- Helper function to predict the reaction outcome ---
    def predict_rearrangement(c1_substituents, c2_substituents):
        """
        Predicts the migrating group based on stability and aptitude rules.
        Returns the name of the migrating group.
        """
        # Calculate stability score for carbocations at C1 and C2
        stab_c1 = sum(groups[g]['stab_score'] for g in c1_substituents)
        stab_c2 = sum(groups[g]['stab_score'] for g in c2_substituents)

        if stab_c1 > stab_c2:
            # Cation forms at C1, so a group from C2 migrates
            migrating_carbon_groups = c2_substituents
        elif stab_c2 > stab_c1:
            # Cation forms at C2, so a group from C1 migrates
            migrating_carbon_groups = c1_substituents
        else:
            # This case is not encountered in this problem
            return "Ambiguous Stability"

        # Determine the migrating group based on highest migratory aptitude
        migrating_group = max(migrating_carbon_groups, key=lambda g: groups[g]['mig_apt'])
        return groups[migrating_group]['name']

    # --- Verify each reaction against the rules ---

    # Reaction A: 3-methyl-4-phenylhexane-3,4-diol
    # Structure: Et-C(OH)(Me) - C(OH)(Ph)-Et
    # C1 substituents: ['Et', 'Me'], C2 substituents: ['Ph', 'Et']
    migrating_group_A = predict_rearrangement(['Et', 'Me'], ['Ph', 'Et'])
    if migrating_group_A != 'ethyl':
        return f"Reasoning for A is wrong: The migrating group should be ethyl, but was predicted as {migrating_group_A}."
    # The predicted product (from ethyl migration) is indeed 3-ethyl-3-phenylpentan-2-one.
    if answer_to_check['A'] != '3-ethyl-3-phenylpentan-2-one':
        return f"Answer for A is wrong: The correct product is '3-ethyl-3-phenylpentan-2-one', but the answer gives '{answer_to_check['A']}'."

    # Reaction B: 3-(4-hydroxyphenyl)-2-phenylpentane-2,3-diol
    # Structure: Me-C(OH)(Ph) - C(OH)(ArOH)-Et
    # C1 substituents: ['Me', 'Ph'], C2 substituents: ['ArOH', 'Et']
    migrating_group_B = predict_rearrangement(['Me', 'Ph'], ['ArOH', 'Et'])
    if migrating_group_B != 'phenyl':
        return f"Reasoning for B is wrong: The migrating group should be phenyl, but was predicted as {migrating_group_B}."
    # The predicted product (from phenyl migration) is 3-(4-hydroxyphenyl)-3-phenylpentan-2-one.
    if answer_to_check['B'] != '3-(4-hydroxyphenyl)-3-phenylpentan-2-one':
        return f"Answer for B is wrong: The correct product is '3-(4-hydroxyphenyl)-3-phenylpentan-2-one', but the answer gives '{answer_to_check['B']}'."

    # Reaction C: 1,1,2-tris(4-methoxyphenyl)-2-phenylethane-1,2-diol
    # Structure: (ArOMe)2-C(OH) - C(OH)(ArOMe)(Ph)
    # C1 substituents: ['ArOMe', 'ArOMe'], C2 substituents: ['ArOMe', 'Ph']
    migrating_group_C = predict_rearrangement(['ArOMe', 'ArOMe'], ['ArOMe', 'Ph'])
    if migrating_group_C != '4-methoxyphenyl':
        return f"Reasoning for C is wrong: The migrating group should be 4-methoxyphenyl, but was predicted as {migrating_group_C}."
    # The predicted product (from 4-methoxyphenyl migration) is 2,2,2-tris(4-methoxyphenyl)-1-phenylethan-1-one.
    if answer_to_check['C'] != '2,2,2-tris(4-methoxyphenyl)-1-phenylethan-1-one':
        return f"Answer for C is wrong: The correct product is '2,2,2-tris(4-methoxyphenyl)-1-phenylethan-1-one', but the answer gives '{answer_to_check['C']}'."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_pinacol_rearrangement_answer()
print(result)