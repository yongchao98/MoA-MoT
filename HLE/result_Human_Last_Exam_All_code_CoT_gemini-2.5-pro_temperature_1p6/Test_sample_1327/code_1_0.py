# The given sequence
s = [2, 11, 23, 51, 119]

print("Let's find the pattern in the sequence:", s)

# Step 1: Establish the primary relationship s(n) = 2*s(n-1) + C
print("\nThe pattern appears to be s(n) = 2 * s(n-1) + C, where C is a changing number.")
print("Let's determine the sequence of C by applying this rule (starting from the 3rd term):")

c1 = s[2] - 2 * s[1]
print(f"{s[2]} = 2 * {s[1]} + {c1}")

c2 = s[3] - 2 * s[2]
print(f"{s[3]} = 2 * {s[2]} + {c2}")

c3 = s[4] - 2 * s[3]
print(f"{s[4]} = 2 * {s[3]} + {c3}")

# Step 2: Analyze the sequence C
c_sequence = [c1, c2, c3]
print("\nThis gives us a new sequence for C:", c_sequence)
print("Now, let's find the pattern in this new sequence by looking at the differences between its terms.")

d1 = c_sequence[1] - c_sequence[0]
print(f"The first difference is {c_sequence[1]} - {c_sequence[0]} = {d1}")

d2 = c_sequence[2] - c_sequence[1]
print(f"The second difference is {c_sequence[2]} - {c_sequence[1]} = {d2}")

# Step 3: Determine the pattern in the differences and predict the next C
print("\nThe differences [4, 12] form a geometric progression with a ratio of 3 (12 / 4 = 3).")
next_d = d2 * 3
print(f"The next difference in this progression is {d2} * 3 = {next_d}.")

next_c = c_sequence[2] + next_d
print(f"Therefore, the next term in the C sequence is {c_sequence[2]} + {next_d} = {next_c}.")

# Step 4: Calculate the final missing term in the original sequence
last_term = s[4]
result = 2 * last_term + next_c
print(f"\nFinally, we can find the missing term in the original sequence:")
print(f"Next Term = 2 * {last_term} + {next_c} = {result}")
