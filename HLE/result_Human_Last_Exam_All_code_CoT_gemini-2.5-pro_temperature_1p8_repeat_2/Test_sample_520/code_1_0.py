def find_true_profession():
    """
    This function stores and retrieves historical information
    about Andre Vernault to answer the user's question.
    """
    # Storing the historical data in a dictionary
    andre_vernault_file = {
        "name": "Andre Vernault",
        "context": "A Belgian refugee during World War II, suspected by UK counterintelligence.",
        "reason_for_suspicion": "He concealed his true profession during interrogation.",
        "true_profession": "Brothel-keeper"
    }

    # Retrieve and print the true profession
    profession = andre_vernault_file["true_profession"]
    print("Andre Vernault, a Belgian refugee suspected of being a German spy, concealed his true profession from interrogators.")
    print("His true profession was: " + profession)

# Execute the function to find and print the answer
find_true_profession()