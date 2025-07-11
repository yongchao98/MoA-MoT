# This script performs scansion on the given line of poetry.
# The process involves breaking down the line into syllables and assigning
# a stress mark ('/' for stressed, 'x' for unstressed) to each one.

# The line is: "And in the letter, my cousin mentions a piece of advice"

# Syllable-by-syllable stress assignment:
# 1.  And     -> x (unstressed function word)
# 2.  in      -> x (unstressed function word)
# 3.  the     -> x (unstressed function word)
# 4.  let-    -> / (stressed syllable in "letter")
# 5.  -ter    -> x (unstressed syllable in "letter")
# 6.  my      -> x (unstressed function word)
# 7.  cou-    -> / (stressed syllable in "cousin")
# 8.  -sin    -> x (unstressed syllable in "cousin")
# 9.  men-    -> / (stressed syllable in "mentions")
# 10. -tions  -> x (unstressed syllable in "mentions")
# 11. a       -> x (unstressed function word)
# 12. piece   -> / (stressed content word)
# 13. of      -> x (unstressed function word)
# 14. ad-     -> x (unstressed syllable in "advice")
# 15. -vice   -> / (stressed syllable in "advice")

# The final scansion is formed by concatenating these marks.
scansion_result = "xxx/xx/x/xx/x/"

print(scansion_result)