import textwrap

def find_hypothetical_word():
    """
    This function explains and determines the Modern English word for 'sister'
    if the Norse influence on the English language had not occurred.
    """

    # 1. The current word "sister" comes from Old Norse "systir".
    # This word replaced the native Old English word due to Viking invasions.
    old_norse_word = "systir"

    # 2. The original Old English (Anglo-Saxon) word was "sweostor".
    # This is a native Germanic word.
    old_english_word = "sweostor"

    # 3. To find the modern equivalent, we trace the likely evolution of "sweostor".
    # - The initial vowel sound 'eo' would likely have simplified to 'e', similar
    #   to how Old English 'seofon' became 'seven'.
    # - The unstressed ending '-or' would likely have been dropped or reduced.
    # This process would result in the Middle English form 'swester'.
    # - This word would likely have survived into Modern English, much like its
    #   German cognate "Schwester" did.
    hypothetical_modern_english_word = "swester"

    # Print the explanation
    explanation = f"""
The modern English word 'sister' is derived from the Old Norse word '{old_norse_word}'. This word was adopted into the language during the Viking Age and eventually replaced the native Old English term.

The original Old English (Anglo-Saxon) word for sister was '{old_english_word}'.

If the language had evolved without this Norse influence, '{old_english_word}' would have followed the natural sound changes of English. It would have likely simplified to a form like '{hypothetical_modern_english_word}', similar to how its cognate in German became 'Schwester'.

Therefore, the Modern English word for 'sister' would most likely be:
    """
    print(textwrap.dedent(explanation).strip())
    print(f"\n{hypothetical_modern_english_word}")

find_hypothetical_word()