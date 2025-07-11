# Step 1: Analyze the truth value of each statement based on logic.
s1_is_logically_true = True  # A tautology because "if [impossible event]..." is always true.
s5_is_logically_false = True # A contradiction because 18 is not older than 19.

# Step 2: Use the "Only one statement is true" rule to determine the truth values of S3 and S4.
# We know S1 is True and S5 is False.
# If S4 were True, S3 would also be True (since equal numbers of males and females makes an even total).
# That would mean S1, S3, and S4 are all True, which violates the "only one is true" rule.
# Therefore, S4 must be False.
s4_is_true = False

# Now, could S3 be True? If S3 were True, then we'd have two True statements (S1 and S3).
# This also violates the "only one is true" rule.
# Therefore, S3 must also be False.
s3_is_true = False

# Step 3: Print the conclusion about which statement is the true one.
print("Step-by-step analysis:")
print(f"Statement S1 is logically True.")
print(f"Statement S5 is logically False.")
print(f"If S4 were true, S3 would also be true, violating the 'only one is true' rule. So, S4 must be False.")
print(f"If S3 were true, there would be two true statements (S1 and S3), violating the rule. So, S3 must be False.")
print("\nConclusion of Analysis:")
print(f"The only true statement is S1.")

# Step 4: Identify the day from the content of the true statement.
# The puzzle states that Fred said the one true statement on the day he tells the truth.
# The true statement is S1: "My name is Fred; and if yesterday was after tomorrow, it would be Wednesday."
# The only clue within this statement that points to a specific day is the mention of "Wednesday".
day_mentioned_in_s1 = "Wednesday"

print(f"Fred said statement S1 on his truth-telling day.")
print(f"The statement is: 'My name is Fred; and if yesterday was after tomorrow, it would be {day_mentioned_in_s1}'.")
print(f"The only piece of information pointing to a specific day in the true statement is the day itself.")
print("\nFinal Answer:")
print(f"The day Fred tells the truth is {day_mentioned_in_s1}.")
