def kalabari_attire_info():
    """
    This function explains the cultural association of the 'etibo'
    dress in Kalabari culture and identifies the correct group.
    """
    # Define the provided answer choices
    choices = {
        'A': 'No special category',
        'B': 'Amanyanabo (sovereign)',
        'C': 'Alabo (chiefs)',
        'D': 'Asawo (young men)',
        'E': 'Opu asawo (gentlemen)'
    }
    
    # The correct answer key based on cultural research
    correct_key = 'C'
    
    # Explanation
    explanation = (
        "The etibo is a formal shirt in Kalabari culture, and its use is "
        "linked to social status. Research into Kalabari traditions reveals "
        "that this particular style of knee-length shirt is most prominently "
        "and formally associated with the chieftaincy class."
    )
    
    # Print the findings
    print(explanation)
    print("\nThe etibo is primarily associated with the following category:")
    print(f"Answer: {correct_key}. {choices[correct_key]}")

# Execute the function to display the answer
kalabari_attire_info()