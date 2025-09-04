import collections

def solve_chemistry_problem():
    """
    This function encapsulates the logic to solve the ring-closing metathesis problem
    and includes tests to verify the solution, following a TDD approach.
    """

    def construct_name_from_subs(subs_dict, base_name):
        """Helper function to build an IUPAC name from substituent info."""
        sub_parts = []
        # Sort by substituent name for alphabetical order in the final name
        for name in sorted(subs_dict.keys()):
            locants = sorted(subs_dict[name])
            prefix = ""
            if len(locants) == 2: prefix = "di"
            if len(locants) == 3: prefix = "tri"
            sub_parts.append(f"{','.join(map(str, locants))}-{prefix}{name}")
        return f"{'-'.join(sub_parts)}-{base_name}"

    def get_rcm_product_name(start_molecule_info):
        """
        Predicts the IUPAC name of the major product from a ring-closing metathesis reaction.
        """
        # 1. Identify inner carbons and determine ring size
        db1_end = start_molecule_info['double_bonds'][0][1]
        db2_start = start_molecule_info['double_bonds'][1][0]
        ring_size = db2_start - db1_end + 1
        ring_base_name = "cyclopent-1-ene" if ring_size == 5 else "cyclohex-1-ene"
        ring_substituents = {p: n for p, n in start_molecule_info['substituents'].items() if db1_end <= p <= db2_start}

        # 2. Evaluate the two possible IUPAC numbering schemes for the product ring
        # Possibility A: New C1=old C6, New C2=old C2
        map_A = {1: db2_start, 2: db1_end}
        for i in range(ring_size - 2): map_A[i + 3] = db1_end + 1 + i
        subs_A = collections.defaultdict(list)
        for new_p, old_p in map_A.items():
            if old_p in ring_substituents: subs_A[ring_substituents[old_p]].append(new_p)
        
        # Possibility B: New C1=old C2, New C2=old C6
        map_B = {1: db1_end, 2: db2_start}
        for i in range(ring_size - 2): map_B[i + 3] = db2_start - 1 - i
        subs_B = collections.defaultdict(list)
        for new_p, old_p in map_B.items():
            if old_p in ring_substituents: subs_B[ring_substituents[old_p]].append(new_p)

        # 3. Apply IUPAC rules to choose the correct numbering
        locants_A = sorted([l for ls in subs_A.values() for l in ls])
        locants_B = sorted([l for ls in subs_B.values() for l in ls])
        
        chosen_subs = None
        if locants_A < locants_B:
            chosen_subs = subs_A
        elif locants_B < locants_A:
            chosen_subs = subs_B
        else:  # Tie-break with alphabetical order
            sorted_names = sorted(subs_A.keys())
            for name in sorted_names:
                if min(subs_A.get(name, [99])) < min(subs_B.get(name, [99])):
                    chosen_subs = subs_A; break
                elif min(subs_B.get(name, [99])) < min(subs_A.get(name, [99])):
                    chosen_subs = subs_B; break
            if not chosen_subs: chosen_subs = subs_A # Default if all else is equal
        
        return construct_name_from_subs(chosen_subs, ring_base_name)

    # Define test cases
    tests = [
        {
            "case_name": "Question's starting material",
            "input": {
                'name': '5-isopropyl-3,4-dimethylocta-1,6-diene',
                'double_bonds': [(1, 2), (6, 7)],
                'substituents': {3: 'methyl', 4: 'methyl', 5: 'isopropyl'}
            },
            "expected_product": "3-isopropyl-4,5-dimethylcyclopent-1-ene",
            "expected_option": "B"
        },
        {
            "case_name": "Test case for a cyclohexene product",
            "input": {
                'name': '4-isopropyl-5,6-dimethylocta-1,7-diene',
                'double_bonds': [(1, 2), (7, 8)],
                'substituents': {4: 'isopropyl', 5: 'methyl', 6: 'methyl'}
            },
            "expected_product": "5-isopropyl-3,4-dimethylcyclohex-1-ene",
            "expected_option": "D"
        }
    ]

    # Run tests
    main_question_test = tests[0]
    try:
        actual_product = get_rcm_product_name(main_question_test['input'])
        assert actual_product == main_question_test["expected_product"]
        
        # Run the second test for robustness
        second_test = tests[1]
        second_product = get_rcm_product_name(second_test['input'])
        assert second_product == second_test['expected_product']

        print(f"All tests passed\n{main_question_test['expected_option']}")

    except AssertionError:
        print(f"Test failed: {main_question_test['case_name']}")
        print(f"Input: {main_question_test['input']['name']}")
        print(f"Expected: {main_question_test['expected_product']}")
        print(f"Got: {actual_product}")
    except Exception as e:
        print(f"An error occurred during testing: {e}")

solve_chemistry_problem()