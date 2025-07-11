import sys
import io

# Set stdout to utf-8 to ensure Arabic characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def solve_manuscript_task():
    """
    This function provides the solution to the manuscript analysis task.
    The analysis of the verse and matres lectionis has been performed beforehand
    based on the provided image and instructions.
    """

    # Task 1: Identify the Verse.
    # The text spans Exodus 7:8 and the beginning of Exodus 7:9.
    # As per the required format (single verse), the starting verse is provided.
    # Format: "Book. Chapter:Verse"
    verse_identification = "Exo. 7:8"

    # Task 2: Compare Matres Lectionis.
    # The comparison results are listed in the order they appear in the Hebrew text (right to left),
    # separated by single spaces.
    # Breakdown of comparison points:
    # 1. וַיֹּאמֶר (BHS: ויאמר) -> ويومر: Substitution of 'aleph' with 'waw' -> או
    # 2. מֹשֶׁה (BHS: משה) -> موسا: Addition of 'waw' for /o/, substitution of 'heh' with 'alif' -> ו הא
    # 3. וְאֶל (BHS: ואל) -> والى: Addition of 'ya' for /e/ -> ي
    # 4. אַהֲרֹן (BHS: אהרן) -> اهرون: Addition of 'waw' for /o/ -> و
    # 5. לֵאמֹר (BHS: לאמר) -> لامور: Addition of 'waw' for /o/ -> و
    # 6. כִּי (BHS: כי) -> كى: Substitution of 'yod' with 'alif maqsura' -> יى
    # 7. פַּרְעֹה (BHS: פרעה) -> فرعو: Substitution of 'heh' with 'waw' -> הו
    # 8. לֵאמֹר (BHS: לאמר) -> لالمور: Addition of 'waw' for /o/ -> و
    matres_comparison = "או ו הא ي و و יى הו ו"

    # Combine the two answers with a single comma and no space as requested.
    final_answer = f"{verse_identification},{matres_comparison}"

    print(final_answer)

solve_manuscript_task()