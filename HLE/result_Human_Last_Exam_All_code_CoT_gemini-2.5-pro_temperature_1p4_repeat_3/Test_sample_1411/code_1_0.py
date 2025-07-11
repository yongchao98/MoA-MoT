def get_watch_inscription():
    """
    This function provides the text written above the date window on the
    Heuer Monaco watch associated with Steve McQueen.

    The watch is the Heuer Monaco, reference 1133.
    The date window is at the 6 o'clock position.
    The text located directly above it is "Automatic Chronograph".
    The function will print this text in all lower case as requested.
    """
    inscription = "Automatic Chronograph"
    
    # Convert the inscription to lower case for the final answer.
    answer = inscription.lower()
    
    print(answer)

get_watch_inscription()