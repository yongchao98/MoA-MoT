# This script prints the solution to the two-part Hebrew phonology question.
# Part 1 identifies the primary stressed syllables in the first seven words of Psalm 74:1.
# Part 2 identifies the last syllable with marked secondary stress in 1 Chronicles 5:10.
# The results are combined into a single string as per the specified format.

syllables_part1 = "כִּיל אָ לָ לֹ נַח נֶ שַׁן"
syllable_part2 = "עָ"

# The final output combines the two parts with a comma and no space.
# It consists only of Hebrew characters, niqqud, te'amim, spaces, and a single comma.
print(f"{syllables_part1},{syllable_part2}")