def provide_answer():
    """
    This function formats the final answer as requested by the user.
    """
    answer_a = "False"
    answer_b = "True"
    answer_c = ["(HA)", "(AS)", "(PA)"]
    
    # Sort the types alphabetically for a consistent output format.
    answer_c.sort()
    
    # Print the final answer in the required format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) [{', '.join(answer_c)}]")

provide_answer()
print("<<<")
provide_answer()
print(">>>")
