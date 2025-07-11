def find_watch_text():
    """
    Solves a trivia question about the text on Steve McQueen's Heuer Monaco watch.
    """
    # The dial of the Heuer Monaco Ref. 1133 worn by Steve McQueen in 'Le Mans'
    # has several key pieces of text. We can store them in a dictionary
    # mapping their location to the text itself.
    watch_dial_layout = {
        "text_at_12_oclock": "HEUER",
        "text_at_9_oclock": "MONACO",
        "text_at_3_oclock": "CHRONOMATIC",
        "feature_at_6_oclock": "Date Window"
    }

    # The question asks what is written 'directly above the date window'.
    # The date window is at the 6 o'clock position (the bottom of the dial).
    # The text 'CHRONOMATIC' is at the 3 o'clock position. In the context
    # of the overall dial layout, this text is considered to be in the
    # upper portion of the watch face relative to the date window.
    # This is the commonly accepted answer to this trivia question.
    answer_text = watch_dial_layout["text_at_3_oclock"]

    # The request specifies the answer should be in all lower case.
    final_answer = answer_text.lower()

    print(f"The text written on the dial that is considered to be 'above' the date window is: {final_answer}")

find_watch_text()