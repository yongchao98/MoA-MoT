def check_name_and_year():
    """
    Addresses the user's question about when Goodluck Jonathan identified as 'Azikiwe'.
    """
    # Goodluck Jonathan's official full name.
    full_name = "Goodluck Ebele Jonathan"

    # The name in question.
    other_name = "Azikiwe"

    # Explanation based on public records.
    explanation = (
        f"Based on public records, {full_name} has never publicly identified himself as '{other_name}'.\n"
        f"His middle name is 'Ebele', not 'Azikiwe'.\n"
        "Therefore, there is no year in which this event occurred."
    )

    print(explanation)

# Run the function to provide the answer.
check_name_and_year()