# The final sequence of numbers is determined by solving the riddle.
# 1. "the fifth... lastly follows" -> The last number is 5.  (?,?,?,?,5)
# 2. "three... protects the last" -> 3 is in the fourth position. (?,?,?,3,5)
# 3. "Number 1 is the first" -> 1 is in the first position. (1,?,?,3,5)
# 4. "Number 4... goes behind [two]" -> The remaining spots are filled by 2 and 4. (1,2,4,3,5)
final_sequence = [1, 2, 4, 3, 5]

# The prompt asks to "output each number in the final equation!".
# I will create a string that represents the sum of the numbers in the sequence.
total = sum(final_sequence)
equation_string = " + ".join(map(str, final_sequence)) + f" = {total}"

print(equation_string)