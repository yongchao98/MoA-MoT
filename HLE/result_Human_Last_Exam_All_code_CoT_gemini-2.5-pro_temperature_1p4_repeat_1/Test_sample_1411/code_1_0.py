def find_watch_text():
    """
    This function stores and retrieves information about the text on the
    face of the Heuer Monaco watch famously worn by Steve McQueen.
    """
    
    # Information about the watch based on known facts and images.
    watch_brand = "Heuer"
    watch_model = "Monaco"
    date_window_position = "6 o'clock"
    
    # Text found on the watch face.
    text_above_date_window = "automatic"
    text_below_main_logo = "monaco"
    
    # The question asks for the text directly above the date window.
    answer = text_above_date_window
    
    print(f"Watch Brand: {watch_brand}")
    print(f"Watch Model: {watch_model}")
    print(f"Location of Date Window: {date_window_position}")
    print(f"Text directly above the date window: {answer}")

find_watch_text()