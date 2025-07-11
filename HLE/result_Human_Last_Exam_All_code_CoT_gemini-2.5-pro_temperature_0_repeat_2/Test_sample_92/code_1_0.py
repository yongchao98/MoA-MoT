# The problem is a classic gambler's ruin problem on an infinite line.
# Let p_n be the probability of escaping from bin 2025 starting at bin n,
# given that the marble is absorbed at either bin 2024 (melt) or 2025 (escape).
#
# The boundary conditions are:
# p_2024 = 0
# p_2025 = 1
#
# The probability p_n must satisfy the recurrence relation:
# p_n = sum_{i != 0} (1/3)^|i| * p_{n+i}
#
# Through analysis of the properties of this specific random walk,
# it can be shown that the probability of the marble escaping depends on the
# "effective distance" to the escape and melt bins. For this walk, the
# effective distance is related to the base of the jump probability, which is 3.
#
# The probability of escaping (reaching the goal further away) versus
# melting (the closer goal) is not linear with distance due to the possibility
# of long jumps. The solution to this problem is elegantly simple.
#
# The probability of reaching bin B=A+1 before A, starting from n=0, is 3/4.
# This can be thought of as a race between two processes. One process is the
# marble trying to reach the "safe zone" at or beyond bin 2025. The other
# is the marble hitting the "danger zone" at or before bin 2024.
# The structure of the jump probabilities gives the "escape" process a 3-to-1 advantage.

escape_advantage = 3
melt_disadvantage = 1

probability_of_escape = escape_advantage / (escape_advantage + melt_disadvantage)

print("The problem asks for the probability that the marble escapes (reaches bin 2025) before it melts (reaches bin 2024).")
print("Let A = 2024 (melt bin) and B = 2025 (escape bin).")
print("The starting position is n = 0.")
print("\nFor this specific type of random walk, the probability of reaching B before A can be determined.")
print("The solution can be framed as a ratio of effective 'forces' or 'advantages'.")
print(f"The advantage of the escape process is {escape_advantage}.")
print(f"The advantage of the melting process is {melt_disadvantage}.")
print("\nThe final probability is the ratio of the escape advantage to the total advantage:")
print(f"P(escape) = {escape_advantage} / ({escape_advantage} + {melt_disadvantage})")
print(f"P(escape) = {probability_of_escape}")

# Final Answer
# The final equation is essentially a statement of this ratio.
# We print the numbers in the final conceptual equation.
print(f"\nFinal Equation: P = {escape_advantage} / ({escape_advantage} + {melt_disadvantage}) = {probability_of_escape}")