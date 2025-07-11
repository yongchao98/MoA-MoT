import sys
import io

# Ensure the output is printed in UTF-8
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def generate_solution():
    """
    This function provides the solution based on the analysis of the provided Judaeo-Arabic text.
    It identifies the author and determines the primary stressed syllables of the first ten words of the text body.
    """
    # 1) The name of the author of the book (one word, in English).
    # The text is Chapter 75 of Maimonides' "Guide for the Perplexed", discussing Kalam proofs for God's unity.
    author_name = "maimonides"

    # 2) A list of the syllables with primary word stress for the first 10 words of the main body of the text.
    # The first 10 words in Arabic are:
    # قَالُوا (qā-lū) -> Stressed: قَا (qā)
    # أَنَّ (ʾan-na) -> Stressed: أَنْ (ʾan)
    # هٰذَا (hā-dhā) -> Stressed: هَا (hā)
    # ٱلَّذِي (al-la-dhī) -> Stressed: أَلْ (ʾal)
    # ذِي (dhī) -> Stressed: ذِي (dhī)
    # ٱلْوُجُودِ (al-wu-jū-di) -> Stressed: جُو (jū)
    # عَلَىٰ (ʿa-lā) -> Stressed: لَا (lā)
    # كَوْنِهِ (kaw-ni-hi) -> Stressed: كَوْ (kaw)
    # صَانِعُهُ (ṣā-ni-ʿu-hu) -> Stressed: صَا (ṣā)
    # وَمُوجِدُهُ (wa-mū-ji-du-hu) -> Stressed: مُو (mū)

    # The stressed syllables in unvocalized Arabic script:
    stressed_syllables = [
        "قا",  # qā
        "ان",  # ʾan
        "ها",  # hā
        "ال",  # ʾal
        "ذي",  # dhī
        "جو",  # jū
        "لا",  # lā
        "كو",  # kaw
        "صا",  # ṣā
        "مو"   # mū
    ]
    
    # Format the final output string
    syllables_string = " ".join(stressed_syllables)
    final_answer = f"{author_name}\n\n{syllables_string}"

    print(final_answer)

generate_solution()