# The user wants to know the 32-character string that should replace the '@' symbol.
# The complete, working one-liner would be:
# from re import*;print(*(sorted({int(id)for id in findall(r'\d+', input())}))or ["NO"])
#
# The part that replaces the '@' placeholder is from the beginning of the print
# arguments up to the start of the findall arguments.
# Let's reconstruct that string.
#
# It starts with the splat operator '*' to handle printing.
# Then `sorted()` to sort the numbers.
# Then a set comprehension `{}` to get unique numbers.
# Inside, `int(id)` converts to integer, and `for id in findall` iterates.
# A 2-character variable 'id' is used instead of 'n' to make the total length 32.
#
# The string is: *(sorted({int(id)for id in findall(
#
# Let's verify the character count:
# *             (1)
# (             (1)
# sorted        (6)
# (             (1)
# {             (1)
# int           (3)
# (             (1)
# id            (2)
# )             (1)
# for           (3)
# id            (2)
# in            (2)
# findall       (7)
# (             (1)
# TOTAL: 1+1+6+1+1+3+1+2+1+3+2+2+7+1 = 32 characters.

print("*(sorted({int(id)for id in findall(")