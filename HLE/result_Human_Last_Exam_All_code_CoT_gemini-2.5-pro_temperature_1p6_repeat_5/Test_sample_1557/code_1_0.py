import json

def evaluate_wellington_statements():
    """
    Analyzes historical statements about the Duke of Wellington's career
    to identify the correct ones.
    """
    # A list of tuples, each containing the statement number and a boolean for its historical accuracy.
    # The reasoning for each True/False evaluation is provided in comments.
    statements_analysis = [
        (1, True),   # Correct. Indian commissariat experience was adapted for the Peninsular War with great success.
        (2, False),  # Incorrect. His intelligence system was famously successful in the Peninsular War (Europe).
        (3, False),  # Incorrect. The British Army was slow to reform; this was not standardized in 1815.
        (4, False),  # Incorrect. Sandhurst was founded much earlier, around 1801-1802, not 1829.
        (5, False),  # Incorrect. Logistical systems and principles were indeed transferable and influential.
        (6, True),   # Correct. Using local auxiliary forces was a model successfully applied in other colonial contexts.
        (7, False),  # Incorrect. This directly contradicts the consensus that his Indian experience was vital for his Peninsular logistics.
        (8, True),   # Correct. The use of 'flying columns' is a clear tactical thread from India to the Peninsula and other colonial wars like the one in Burma.
        (9, False),  # Incorrect. The 1813 Charter Act was driven by broad political/economic factors, not directly by Wellington's military principles.
        (10, False)  # Incorrect. His Indian experience in civil administration was highly relevant and applied in the Peninsula.
    ]

    correct_options = []
    for number, is_correct in statements_analysis:
        if is_correct:
            correct_options.append(number)

    # Sort the numbers in ascending order
    correct_options.sort()

    # Format the result as a comma-separated string
    result_string = ", ".join(map(str, correct_options))

    # Print the final result as requested. Each number is output in the final line.
    print(result_string)

if __name__ == "__main__":
    evaluate_wellington_statements()