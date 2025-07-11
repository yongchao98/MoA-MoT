import re

def final_answer():
    """
    This function prints the answers to the user's questions,
    and also prints the symbolic and numeric components of the resulting expressions.
    """
    # The derived answers to the questions
    ans_a = "w^d = 0"
    ans_b = "w^d * r"
    ans_c = "yes"

    # Print the final answer in the specified format
    print(f"(a) [{ans_a}] (b) [{ans_b}] (c) [{ans_c}]")

    # Fulfilling the requirement to output each number in the final equation.
    # For (a), the number is 0. The other elements are symbolic variables.
    numbers_a = re.findall(r'\d+', ans_a)
    print("\nIn equation (a), the number is:", numbers_a[0])

final_answer()