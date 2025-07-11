def find_watch_inscription():
    """
    This function stores and retrieves information about the watch in question
    to answer the user's query.
    """
    # Details of the watch based on public information about the auction
    watch = {
        "brand": "Heuer",
        "model": "Monaco",
        "reference": "1133",
        "association": "Steve McQueen in 'Le Mans'",
        "feature_location": "Directly above the date window at 6 o'clock",
        "inscription": "chronograph"
    }

    # Retrieve the inscription and convert it to lower case as requested
    answer = watch["inscription"].lower()

    print(f"On the Heuer Monaco watch, the text written directly above the date window is: {answer}")

find_watch_inscription()