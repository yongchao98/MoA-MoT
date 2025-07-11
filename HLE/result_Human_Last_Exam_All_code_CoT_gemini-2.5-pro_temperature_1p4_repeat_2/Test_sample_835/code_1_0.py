import sys
import io

# Ensure the output is encoded in UTF-8 to correctly display Hebrew characters.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Part 1: A list of the syllables carrying primary word stress in the first seven words of the verse.
# The syllables are identified based on Tiberian Hebrew accentuation rules.
stressed_syllables = "כִּ֗יל סָ֥ף מָ֣ה לֹ֭ נַ֣ח נֶ֑ שַׁ֥ן"

# Part 2: The last occurrence of a syllable with secondary stress (marked by a meteg) in the second verse.
last_secondary_stress_syllable = "עָֽ"

# The final answer combines both parts, separated by a comma with no space,
# as per the instructions.
final_answer = f"{stressed_syllables},{last_secondary_stress_syllable}"

print(final_answer)