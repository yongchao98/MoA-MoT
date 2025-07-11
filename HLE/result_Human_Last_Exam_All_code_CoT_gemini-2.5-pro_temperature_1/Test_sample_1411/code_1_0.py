def find_watch_text():
    """
    This function identifies and prints the text found above the date window
    on the watch famously worn by Steve McQueen.

    The watch is the Heuer Monaco (ref. 1133). On this model, the date window
    is located at the 6 o'clock position. The text directly above the window
    is the model name, "Monaco".
    """
    
    # The text on the watch face.
    text_on_watch = "Monaco"
    
    # Convert the text to the requested all-lowercase format.
    answer = text_on_watch.lower()
    
    # Print the final answer.
    print(answer)

find_watch_text()