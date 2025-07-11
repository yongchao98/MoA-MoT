def trace_etymology():
    """
    This function stores and prints the etymological development of a PIE root
    through several Germanic languages.
    """
    pie_form = "*seh2gieti"
    meaning = "'(s)he's giving a sign'"

    # Derived forms based on historical linguistic rules
    proto_germanic = "*sōkijaną"  # from PIE *sāg-eyo-
    old_norse = "sœkja"
    old_english = "sēcan"
    old_high_german = "suohhan"

    # Print the final result in a clear sentence
    # The prompt asked to output each word, which this print statement does.
    print(
        f"The Proto-Indo-European form {pie_form} {meaning} gives us the following cognates:\n"
        f"Proto-Germanic: {proto_germanic}\n"
        f"Old Norse: {old_norse}\n"
        f"Old English: {old_english}\n"
        f"Old High German: {old_high_german}"
    )

trace_etymology()