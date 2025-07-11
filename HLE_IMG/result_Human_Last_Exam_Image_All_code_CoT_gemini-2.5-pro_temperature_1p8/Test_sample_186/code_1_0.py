import sys
import io

# Ensure the output is printed in UTF-8 to correctly handle Arabic and Hebrew characters.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def generate_manuscript_analysis():
    """
    This function contains the analysis of the manuscript text as per the user's request.
    The analysis is hardcoded based on the detailed examination of the manuscript image and
    comparison with the Biblia Hebraica Stuttgartensia (BHS).

    Verse Identification: The text corresponds to Exodus 11:9.
    BHS Text: וַיֹּאמֶר יְהוָה אֶל-מֹשֶׁה לֹא-יִשְׁמַע אֲלֵיכֶם פַּרְעֹה לְמַעַן רְבוֹת מוֹפְתַי בְּאֶרֶץ מִצְרָיִם
    Manuscript Transcription (lines 5-6): ويومر ربنا ال موسا لا يسمع اليكم فرعو لمعن ربوا موفتاي في ارض مصرام

    Comparison of Matres Lectionis (right to left):
    1. וַיֹּאמֶר (vayo'mer) -> ويامر (wayāmir): BHS has defective spelling for the 'a' vowel. The Arabic transcription adds an alif (ا) as a mater lectionis. -> Add: ا
    2. מֹשֶׁה (Mosheh) -> موسا (Mūsā): The final Hebrew mater ה is substituted with the Arabic mater ا. -> Substitute: הا
    3. פַּרְעֹה (Par'oh) -> فرعو (Firʿaw): The final Hebrew mater ה is substituted with the Arabic mater و. -> Substitute: הו
    4. רְבוֹת (rəvot) -> ربوا (rubūwā): The transcription represents the feminine plural ending -ot with -ūwā, adding a final mater ا not present in the Hebrew. -> Add: ا
    5. מוֹפְתַי (mofetay) -> موفتاي (mōfatāy): The Hebrew -tay is spelled defectively. The Arabic transcription adds an alif (ا) for the 'a' vowel before the yā'. -> Add: ا
    6. בְּאֶרֶץ (bə'erets) -> في ارض (fī arḍ): The Hebrew prefix bə- is represented by the Arabic word fī, which includes a mater lectionis ي not present in the Hebrew. -> Add: ي
    7. מִצְרָיִם (Mitsrayim) -> مصرام (Miṣrām): The Hebrew mater י in the dual ending -ayim is substituted with the Arabic mater ا in the variant ending -ām. -> Substitute: يا
    """

    verse_identification = "Exo. 11:9"
    matres_comparison = "ا הا הו ا ا ي يا"

    # The final answer combines both parts with a single comma, no space.
    final_answer = f"{verse_identification},{matres_comparison}"
    print(final_answer)

generate_manuscript_analysis()