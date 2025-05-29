from math import comb

# Case 1: All 4 balls in one box
case1 = 1

# Case 2: 3 balls in one box, 1 in the other
case2 = comb(4, 3)

# Case 3: 2 balls in each box
case3 = comb(4, 2) // 2

# Total number of ways
total_ways = case1 + case2 + case3

print(total_ways)