def get_answer():
    """
    This function prints the final answers to the questions.
    """
    answer_a = "\\emptyset, \\mathbb{R}"
    answer_b = "No"
    answer_c = "Yes"
    
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]")

get_answer()
print("<<<")
print(f"(a) {get_answer.__doc__}") # This line is just to provide the required format, not a real execution
print(">>>")