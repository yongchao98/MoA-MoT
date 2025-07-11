import sys
import io

# Ensure the output is printed in UTF-8
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def generate_solution():
    """
    This function generates the final answer based on the analysis of the manuscript.
    The analysis itself is performed based on visual inspection of the manuscript and knowledge
    of Hebrew and Karaite Arabic transcriptions. The code here assembles and prints the result.
    """
    
    # Task 1: Identify the Verse
    # Analysis of lines 5 and 6 on the right-hand page reveals the text is from the book of Numbers.
    # The phrase "wa-yaʿśū lāhīm ṣīṣāh" (and they shall make for themselves fringes) is a clear
    # indicator of Numbers 15:38.
    verse_identification = "Num. 15:38"

    # Task 2: Compare Matres Lectionis
    # Comparing the BHS text of Numbers 15:38 with the manuscript's transcription word-by-word (right-to-left).
    # The BHS text for the analyzed portion is:
    # דַּבֵּר אֶל-בְּנֵי יִשְׂרָאֵל ...וְאָמַרְתָּ אֲלֵהֶם וְעָשׂוּ לָהֶם צִיצִת עַל-כַּנְפֵי בִגְדֵיהֶם לְדֹרֹתָם
    # The MS transcription is:
    # هادور الى بني يسرائل...وقول تيهم وياعسوا لاهيم صيصاة عل كنف بكاديهيم لدوروتهيم

    comparison_items = [
        ("דַּבֵּר", "هادور", "ا ו"),  # Hebrew has no matres; Arabic adds 'alif' and 'waw'.
        ("אֶל", "الى", "ا"),        # Hebrew has no mater; Arabic adds 'alif maqsura'.
        ("בְּנֵי", "بني", ""),       # Match: י -> ي
        ("יִשְׂרָאֵל", "يسرائل", "ا ي"), # Heb 'rā' (no mater), 'ēl' (no mater); Arb 'rā' (ا), 'īl' (ي).
        ("אֲלֵהֶם", "...تيهم", "ي"), # Heb 'ē' (no mater); Arb 'ī' (ي).
        ("וְעָשׂוּ", "وياعسوا", "ا"), # Heb 'ā' (no mater); Arb 'ā' (ا). 'ū' -> 'ū' matches.
        ("לָהֶם", "لاهيم", "ا ي"), # Heb 'ā' (no mater), 'e' (no mater); Arb 'ā' (ا), 'ī' (ي).
        ("צִיצִת", "صيصاة", "ا"),     # Heb 'i' (no mater); Arb 'ā' (ا). 'ī' -> 'ī' matches.
        ("עַל", "عل", ""),          # Match: no matres.
        ("כַּנְפֵי", "كنف", "noי"),   # Hebrew has 'י'; Arabic does not.
        ("בִּגְדֵהֶם", "بكاديهيم", "ا ي"), # Heb 'bi', 'dē' (no matres); Arb 'kā' (ا), 'dī' (ي).
        ("לְדֹרֹתָם", "لدوروتهيم", "و ي") # Heb 'ō', 'ā' (no matres); Arb 'ō' (و), 'ī' (ي).
    ]
    
    # Building the final string from the non-empty comparison results.
    matres_list = []
    for hebrew, arabic, result in comparison_items:
        if result:
            # The split() is to handle cases like "ا ו" correctly.
            matres_list.extend(result.split())

    matres_string = " ".join(matres_list)

    final_answer = f"{verse_identification},{matres_string}"
    
    print(final_answer)

generate_solution()