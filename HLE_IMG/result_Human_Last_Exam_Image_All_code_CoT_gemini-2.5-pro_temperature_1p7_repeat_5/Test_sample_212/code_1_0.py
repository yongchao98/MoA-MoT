def solve():
    """
    This function formats the famous message on the back of Paul Newman's Rolex Daytona.
    """
    # The original message engraved on the watch
    message = "DRIVE CAREFULLY ME"

    # Convert to lowercase
    formatted_message = message.lower()

    # The prompt doesn't ask for removing spaces, so we'll keep them to preserve the message.
    # The message has no punctuation to remove.

    print(formatted_message)

solve()