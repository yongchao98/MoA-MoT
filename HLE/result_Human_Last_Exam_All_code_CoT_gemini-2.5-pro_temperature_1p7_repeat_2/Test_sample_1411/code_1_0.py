def find_watch_text():
    """
    This function identifies and prints the text on a specific watch.
    
    1. The watch is identified as the Heuer Monaco, famously worn by Steve McQueen, which has a date window.
    2. On this watch's dial, the date window is at the 6 o'clock position.
    3. The text "Automatic Chronograph" is located directly above the date window.
    4. The text is converted to lowercase as per the request.
    """
    
    # The text written above the date window on the Heuer Monaco watch.
    text_on_watch = "Automatic Chronograph"
    
    # Convert the answer to all lower case for the final output.
    final_answer = text_on_watch.lower()
    
    print(final_answer)

find_watch_text()