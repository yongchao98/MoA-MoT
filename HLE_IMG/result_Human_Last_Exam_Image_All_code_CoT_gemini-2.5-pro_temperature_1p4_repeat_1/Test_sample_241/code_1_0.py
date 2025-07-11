import sys
import io

# Set stdout to utf-8 to ensure Arabic characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def provide_solution():
    """
    This function prints the solution to the task, including the author's name
    and the stressed syllables from the provided text, as per the user's instructions.
    """
    
    # 1. The name of the author of the book (one word, in English, lowercase).
    author_name = "maimonides"

    # 2. A list of the syllables with primary word stress for the first 10 words.
    # The words are: أَنَا (anā), أُبَيِّنُ (ubayyinu), لَكَ (laka), أَيْضًا (ayḍan),
    # فِي (fī), هَذَا (hādhā), الفَصْلِ (al-faṣli), دَلَائِلَ (dalāʾila),
    # التَّوْحِيدِ (al-tawḥīdi), عَلَى (ʿalā).
    # The stressed syllables are: nā, bay, lak, ʾay, fī, hā, faṣl, lā, ḥīd, lā.
    # The syllables are presented in unvocalized Arabic script, separated by spaces.
    stressed_syllables_list = ["نا", "بي", "لك", "اي", "في", "ها", "فصل", "لا", "حيد", "لا"]
    stressed_syllables_string = " ".join(stressed_syllables_list)
    
    # Print the results in the required format.
    print(author_name)
    print()  # Exactly one line space between sections.
    print(stressed_syllables_string)

provide_solution()

# The final answer is wrapped in <<<>>>
final_answer_string = "maimonides\n\nنا بي لك اي في ها فصل لا حيد لا"
print(f"<<<{final_answer_string}>>>")