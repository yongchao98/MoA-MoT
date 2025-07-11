# Step 1: Define the variables based on the clues.
# X1: Number of deputies in the Russian State Duma.
X1 = 450
# X2: From a pun on the Russian word for barber (парикмахер), which contains "пара" (pair/two).
X2 = 2
# X3: The Ten Commandments.
X3 = 10
# X4: The length in kilometers of the queue in the Bernard Shaw anecdote.
X4 = 8
# X6: The 11th letter 'K' (Potassium, essential for life).
X6 = 11
# X7: From the fairy tale "Seven at one blow".
X7 = 7
# X8: Represents duality.
X8 = 2
# X9: From the Russian proverb "Don't have 100 rubles, have 100 friends".
X9 = 100
# X10: From Kozma Prutkov's aphorism #13 on lying.
X10 = 13
# X11: From the "one Belza" unit of measurement.
X11 = 1

# Step 2: Solve the first equation for X5.
# The equation is: X1*X2 + (X3*X4 - X5*X6) + X7*X4 = X8*X9
# Rearranging to solve for X5:
# X5*X6 = X1*X2 + X3*X4 + X7*X4 - X8*X9
X5 = (X1 * X2 + X3 * X4 + X7 * X4 - X8 * X9) / X6
X5 = int(X5)

# Step 3: Interpret the second line as an equation to find X12.
# The line is: X8's X9's X10 X11 X12.
# This implies a calculation: X12 = X8 * X9 * X10 * X11
X12 = X8 * X9 * X10 * X11

# Step 4: Print the final equation and the result.
# The prompt asks to output each number in the final equation.
print(f"The final calculation to find X12 is:")
print(f"{X8} * {X9} * {X10} * {X11} = {X12}")

# The final answer for X12
# print(f"The value of X12 is {X12}")