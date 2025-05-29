# Possible numbers and letters based on the feedback
possible_numbers = set(range(10))
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Eliminate numbers based on feedback
# 58TD: both numbers incorrect
possible_numbers -= {5, 8}
# 19UI: both numbers incorrect
possible_numbers -= {1, 9}
# 75AC: both numbers incorrect
possible_numbers -= {7, 5}
# 08LP: both numbers incorrect
possible_numbers -= {0, 8}
# 51PB: both numbers incorrect
possible_numbers -= {5, 1}
# 98OJ: both numbers incorrect
possible_numbers -= {9, 8}

# Eliminate letters based on feedback
# 58TD: both letters incorrect
possible_letters -= {'T', 'D'}
# 52SK: both letters incorrect
possible_letters -= {'S', 'K'}
# 19UI: both letters incorrect
possible_letters -= {'U', 'I'}
# 75AC: both letters incorrect and too early
possible_letters -= {'A', 'C'}
# 86WD: both letters incorrect
possible_letters -= {'W', 'D'}
# 72KL: both letters incorrect and too early
possible_letters -= {'K', 'L'}
# 08LP: both letters incorrect
possible_letters -= {'L', 'P'}
# 51PB: both letters incorrect
possible_letters -= {'P', 'B'}
# 96FX: both letters incorrect
possible_letters -= {'F', 'X'}
# 02HU: both letters incorrect
possible_letters -= {'H', 'U'}
# 96FJ: both letters incorrect and too early
possible_letters -= {'F', 'J'}
# 12KC: both letters incorrect and too early
possible_letters -= {'K', 'C'}
# 92TG: both letters incorrect
possible_letters -= {'T', 'G'}
# 98OJ: both letters incorrect
possible_letters -= {'O', 'J'}
# 06ML: both letters incorrect and too early
possible_letters -= {'M', 'L'}
# 02VX: one letter correct but in wrong position; one letter incorrect and too late
possible_letters -= {'V'}

# Remaining possible numbers and letters
print("Possible numbers:", possible_numbers)
print("Possible letters:", possible_letters)