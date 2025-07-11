import textwrap

def explain_hypothetical_word():
    """
    Explains and prints the likely Modern English word for "sister"
    without Norse influence.
    """

    # Step 1: Explain the origin of the current word "sister"
    print("Step 1: The Origin of 'Sister'")
    explanation1 = """\
    The modern English word 'sister' is not a direct descendant from Old English.
    It was borrowed from the Old Norse word 'systir' during the period of
    Viking invasions and settlement in England (roughly 8th to 11th centuries).
    The Norse word replaced the native English one.
    """
    print(textwrap.dedent(explanation1))

    # Step 2: Identify the original Old English word
    print("Step 2: The Native Old English Word")
    explanation2 = """\
    The original word for sister in Old English was 'sweostor'. This word is part
    of the common Germanic heritage, with cognates like 'Schwester' in German
    and 'zuster' in Dutch.
    """
    print(textwrap.dedent(explanation2))

    # Step 3: Trace the hypothetical evolution
    print("Step 3: How 'sweostor' Would Have Evolved")
    explanation3 = """\
    Without Norse influence, 'sweostor' would have undergone natural sound changes:
    - The 'w' sound would likely have been lost, as seen in 'sword' (from 'sweord').
    - The 'eo' vowel combination would have simplified.
    - In many dialects, 'sweostor' became 'suster' in Middle English. This form
      was common and used by writers like Chaucer.
    - The vowel 'u' in 'suster' (pronounced like the 'u' in 'put') would have
      naturally evolved into the modern 'uh' sound (like the 'u' in 'cup').
    """
    print(textwrap.dedent(explanation3))

    # Step 4: State the final conclusion
    print("Step 4: The Resulting Modern English Word")
    final_word = "suster"
    explanation4 = f"""\
    Following this linguistic path, the Middle English 'suster' would have
    become the Modern English word '{final_word}'. It would be pronounced
    to rhyme with 'duster' or 'cluster'.
    """
    print(textwrap.dedent(explanation4))
    print("-----------------------------------------")
    print(f"The hypothetical Modern English word is: {final_word}")


if __name__ == "__main__":
    explain_hypothetical_word()