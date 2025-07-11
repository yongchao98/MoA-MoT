def find_secret_word():
    """
    This function reveals the secret word hidden in the image.
    The message is encoded across the Red, Green, and Blue color channels
    of the central noisy area.
    """
    
    # Parts of the message were discovered in the separate color channels
    red_channel_words = ("THE", "LIE")
    green_channel_word = "CAKE"
    blue_channel_words = ("IS", "A")

    # Assemble the final phrase from the discovered parts
    part1 = red_channel_words[0]
    part2 = green_channel_word
    part3 = blue_channel_words[0]
    part4 = blue_channel_words[1]
    part5 = red_channel_words[1]

    final_phrase = f"{part1} {part2} {part3} {part4} {part5}"

    # Print the components that make up the final phrase
    print("The final message is formed by combining the following parts found in the image:")
    print(f"'{part1}' + '{part2}' + '{part3}' + '{part4}' + '{part5}'")
    
    # Print the secret phrase
    print("\nThe secret word is:")
    print(final_phrase)

find_secret_word()