import sys
import io

# Set stdout to utf-8 to handle Arabic and Hebrew characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def solve_manuscript_puzzle():
  """
  This function provides the solution to the manuscript analysis task.
  It identifies the biblical verse and compares the use of matres lectionis
  between the Hebrew source (BHS) and the Arabic transcription.
  """

  # Task 1: Identify the Verse
  # The text spans Exodus 7:8 and the beginning of Exodus 7:9.
  # Following the requested format, the starting verse is provided.
  verse_identification = "Exo. 7:8"

  # Task 2: Compare Matres Lectionis
  # The findings are listed in order of appearance in the text (right-to-left).
  # 1. וַיֹּאמֶר vs. وايو مر -> Additional و
  # 2. מֹשֶׁה vs. يمو سا -> Additional و, Substitute הא
  # 3. יְדַבֵּר vs. كنار -> Additional ا
  # 4. אֲלֵכֶם vs. العام -> Substitute يا
  # 5. פַּרְעֹה vs. في عو -> Substitute הו
  matres_comparison = "و و הא ا יא הו"

  # Combine the answers into the final format, separated by a single comma.
  final_answer = f"{verse_identification},{matres_comparison}"
  
  print(final_answer)

solve_manuscript_puzzle()