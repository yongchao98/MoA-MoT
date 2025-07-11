# Step 1: Set Y1 based on Clue 1 (Charles Martin Hall's patent year for aluminum process)
Y1_str = "1886"
X1 = int(Y1_str[0])
X2 = int(Y1_str[1])
X3 = int(Y1_str[2])
X4 = int(Y1_str[3])

# Step 2: The numerical string derived from Clue 2
numerical_string = "238438653790"

# Step 3: Parse the numerical string into Y2, Y3, and Y4
Y2_str = numerical_string[0:5]
Y3_str = numerical_string[5:9]
Y4_str = numerical_string[9:12]

# Step 4: Solve for the remaining X variables using the knowns and parsed strings
# Y2 = X5 X6 X2 X7 X6 = X5 X6 8 X7 X6
# From Y2_str = "23843", we can deduce X5, X6, X7
X5 = int(Y2_str[0])  # 2
X6 = int(Y2_str[1])  # 3
X7 = int(Y2_str[3])  # 4
# We check for consistency: X2 is indeed the middle digit, and X6 repeats.
# Y2_str[2] -> 8 (matches X2)
# Y2_str[4] -> 3 (matches X6)

# Y3 = X3 X4 X8 X6 = 8 6 X8 X6
# From Y3_str = "8653", we can deduce X8
X8 = int(Y3_str[2]) # 5
# We check for consistency: X3, X4 are the first two digits, and X6 is the last.
# Y3_str[0] -> 8 (matches X3)
# Y3_str[1] -> 6 (matches X4)
# Y3_str[3] -> 3 (matches X6)

# Y4 = X9 X10 X11
# From Y4_str = "790", we deduce X9, X10, X11
X9 = int(Y4_str[0])
X10 = int(Y4_str[1])
X11 = int(Y4_str[2])

# Step 5: Print the final calculated equations
print(f"Y1 = {X1}{X2}{X3}{X4}")
print(f"Y2 = {X5}{X6}{X2}{X7}{X6}")
print(f"Y3 = {X3}{X4}{X8}{X6}")
print(f"Y4 = {X9}{X10}{X11}")
