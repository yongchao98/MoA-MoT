# Step 1: Decode Clue 1 to establish X1, X2, X3, X4.
# Clue 1 points to Abraham Gesner's 1854 patent for kerosene.
# Y1 = 1854
X1 = 1
X2 = 8
X3 = 5
X4 = 4

# Step 2: Decode Clue 2 to find the remaining X variables.
# The clue suggests taking a famous number sequence (Pi), and removing the digit '9' (pun on 'ge').
# The digits of Pi are: 3.1415926535...
# After removing the '9', the sequence of digits used is: 3, 1, 4, 1, 5, 2, 6, 5, 3, 5...
X5 = 3
X6 = 1
X7 = 4
X8 = 1
X9 = 5
X10 = 2
X11 = 6

# Step 3: Calculate the Y indices based on the decoded X variables.
# The formulas define the Y values as concatenations of the X digits.
Y1_str = f"{X1}{X2}{X3}{X4}"
Y2_str = f"{X5}{X6}{X2}{X7}{X6}"
Y3_str = f"{X3}{X4}{X8}{X6}"
Y4_str = f"{X9}{X10}{X11}"

Y1 = int(Y1_str)
Y2 = int(Y2_str)
Y3 = int(Y3_str)
Y4 = int(Y4_str)

# Step 4: Print the final results in the specified format.
# The output shows each Y value and the equation that forms it.
print(f"Y1 ({Y1}) = X1({X1})X2({X2})X3({X3})X4({X4})")
print(f"Y2 ({Y2}) = X5({X5})X6({X6})X2({X2})X7({X7})X6({X6})")
print(f"Y3 ({Y3}) = X3({X3})X4({X4})X8({X8})X6({X6})")
print(f"Y4 ({Y4}) = X9({X9})X10({X10})X11({X11})")
