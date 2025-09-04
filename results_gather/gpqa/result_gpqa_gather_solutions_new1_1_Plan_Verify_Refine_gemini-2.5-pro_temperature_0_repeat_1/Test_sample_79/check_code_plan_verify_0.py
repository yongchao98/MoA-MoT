def check_answer():
    """
    Checks the logic that the correct answer is the one that provides a unique
    positive integer solution for the values of A, C, G, and T.
    """
    options = {
        'A': 351,
        'B': 333,
        'C': 315,
        'D': 185
    }
    
    final_answer_value = 185 # This corresponds to the provided answer 'D'
    
    solution_counts = {}
    found_solutions = {}

    # Iterate through each multiple-choice option
    for option_letter, option_value in options.items():
        count = 0
        solutions_for_option = []
        
        # Constraints for positive integers:
        # A = 115 - 2*G > 0  => 2*G < 115 => G <= 57
        # C = 61 - 2*T > 0   => 2*T < 61  => T <= 30
        for G in range(1, 58):
            A = 115 - 2 * G
            if A <= 0: continue # A must be positive

            for T in range(1, 31):
                C = 61 - 2 * T
                if C <= 0: continue # C must be positive

                # Calculate the value for the target string 'ACAGTGACC'
                # 3*A + 3*C + 2*G + T
                calculated_value = 3 * A + 3 * C + 2 * G + T
                
                if calculated_value == option_value:
                    count += 1
                    solutions_for_option.append({'A': A, 'C': C, 'G': G, 'T': T})

        solution_counts[option_letter] = count
        found_solutions[option_letter] = solutions_for_option

    # Check the core claim of the provided answer's reasoning
    # The claim is that only option D (185) has a unique solution.
    
    # Verify the final answer value (185)
    if final_answer_value not in options.values():
        return f"The final answer value {final_answer_value} is not one of the options."

    final_answer_letter = [k for k, v in options.items() if v == final_answer_value][0]

    # Check if the final answer's option is indeed unique
    if solution_counts[final_answer_letter] != 1:
        return f"The reasoning is incorrect. The chosen answer {final_answer_letter} ({final_answer_value}) does not have a unique solution. It has {solution_counts[final_answer_letter]} solutions."

    # Check if any other option is also unique
    for letter, count in solution_counts.items():
        if letter != final_answer_letter and count == 1:
            return f"The reasoning is incorrect. Option {letter} also has a unique solution, making the uniqueness argument invalid."

    # Verify the specific unique solution cited in the reasoning
    expected_solution = {'A': 1, 'C': 15, 'G': 57, 'T': 23}
    if found_solutions[final_answer_letter][0] != expected_solution:
        return f"The code found a unique solution for {final_answer_value}, but it is {found_solutions[final_answer_letter][0]}, not the one cited in the reasoning {expected_solution}."

    # If all checks pass, the reasoning and answer are correct.
    return "Correct"

# Run the check
result = check_answer()
if result == "Correct":
    print("Correct")
else:
    print(f"Incorrect. Reason: {result}")
