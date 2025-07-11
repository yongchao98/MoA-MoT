# Lower and upper bounds for 'a' based on n=4
a_low = 3**(1/3)
a_high = 5**(1/4)

# Calculate the range for a^5
a5_low = a_low**5
a5_high = a_high**5

print(f"For a in [{a_low:.4f}, {a_high:.4f}),")
print(f"a^5 must be in [{a5_low:.4f}, {a5_high:.4f}).")

# Find possible integer values for floor(a^5)
first_possible_floor = int(a5_low) + 1
last_possible_floor = int(a5_high)

possible_floors = range(first_possible_floor, last_possible_floor + 1)
print(f"Possible integer values for floor(a^5) are: {list(possible_floors)}.")

# Check if any of these satisfy the n=5 condition (floor(a^5) % 3 == 2)
found_solution = False
for k in possible_floors:
    print(f"Checking k = {k}: {k} % 3 = {k % 3}.")
    if k % 3 == 2:
        found_solution = True

if not found_solution:
    print("\nNone of the possible integer values for floor(a^5) satisfy the condition k % 3 == 2.")
    print("This construction starting with floor(a)=1 fails at n=5.")

# The problem asks for just the final yes/no answer, separated by a comma.
# Part 1 (mod 2): Yes
# Part 2 (mod 3): No
print("\nFinal Answer:")
print("Yes,No")
