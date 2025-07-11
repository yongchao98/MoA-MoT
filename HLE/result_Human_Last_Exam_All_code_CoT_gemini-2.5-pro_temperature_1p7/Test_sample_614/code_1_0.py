# Based on the riddle's clues, we can deduce the order of the five numbers.
# 1. "Number 1 is the first..." -> The sequence starts with 1.
# 2. "...the fifth, who...lastly follows." -> The number 5 is at the end.
# 3. "Number three...protects the last;" -> The number 3 is in the fourth position, before 5.
# 4. "Number 4...always goes behind [Number 2]." -> The numbers 2 and 4 appear as a block (2, 4), filling the remaining second and third spots.
# This gives the sequence 1, 2, 4, 3, 5.

# Let's verify the "bloodline" clue: "the bloodline rules them all but the fifth".
# This suggests a rule for the first four numbers: 1, 2, 4, 3.
# A simple rule is: first number + second number = fourth number.
# In our sequence: 1 + 2 = 3. This rule holds. The fifth number, 5, is not involved.

# Final sequence
num1 = 1
num2 = 2
num3 = 4
num4 = 3
num5 = 5

# The prompt mentions a "final equation", which we interpret as displaying the ordered numbers.
print(f"{num1}, {num2}, {num3}, {num4}, {num5}")
