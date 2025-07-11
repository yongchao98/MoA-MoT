# Step 1: Define the values for X1 through X8, X10, and X11 based on the clues.
X1 = 150  # From Vladimir Mayakovsky's poem "150,000,000".
X2 = 15   # From the aria in "The Barber of Seville", "quindici" (fifteen).
X3 = 500  # From the M500 grade of Portland cement and the 500 brethren in the Bible.
X4 = 32   # From the number of letters in the Russian alphabet in 1931.
X5 = 40   # From the Russian pun on "Magpies" (Сороки) and "forty" (сорок).
X6 = 6    # From the Latin "sex", linking to life and focus.
X7 = 42   # From the duration of the sleeping potion in "Romeo and Juliet".
X8 = 2    # From the two sides (Army and Navy) in the golf match.
X10 = 2   # From Kozma Prutkov's aphorism about two celestial bodies (sun and moon).
X11 = 1   # From Svyatoslav Belza's "bel-etagè", representing the first floor.

# Step 2: Calculate the left side of the first equation to find the product of X8 and X9.
equation1_result = X1 * X2 + (X3 * X4 - X5 * X6) + X7 * X4
X8_X9_product = equation1_result

# Step 3: Calculate X9.
X9 = X8_X9_product / X8

# Step 4: Interpret the second line as an equation and solve for X12.
# X8 * X9 = X10 * X11 * X12
X12 = X8_X9_product / (X10 * X11)

# Step 5: Print the results as requested.
print("Solving the first equation:")
print(f"{X1} * {X2} + ({X3} * {X4} - {X5} * {X6}) + {X7} * {X4} = {X8} * {int(X9)}")
print(f"Result of the left side is: {equation1_result}")
print(f"So, {X8} * X9 = {equation1_result}")
print(f"This gives X9 = {int(X9)}")
print("\nSolving for X12 using the second equation interpretation (X8 * X9 = X10 * X11 * X12):")
print(f"{int(X8_X9_product)} = {X10} * {X11} * X12")
print(f"X12 = {int(X12)}")
