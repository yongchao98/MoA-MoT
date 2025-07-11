# Define the parameters of the problem based on the description.
num_cols = 2023
num_monsters = 2022

# The problem states there is exactly one monster in each row from 1 to 2022,
# and each column contains at most one monster.
# This implies that the 2022 monsters must be in 2022 different columns.
# These 2022 columns are "unsafe".
num_unsafe_cols = num_monsters

# To guarantee reaching the last row, Turbo needs a strategy that works even in the worst-case scenario.
# The worst case is that the single "safe column" (the one with no monster) is the last one Turbo checks.
# A single attempt can, at best, prove that one column is unsafe (by finding a monster in it).
# Therefore, in the worst case, Turbo must use one attempt to find the monster in each of the unsafe columns.
max_failed_attempts = num_unsafe_cols

# After all unsafe columns have been identified through failed attempts,
# the single remaining column is known to be safe. Turbo needs one more attempt
# to travel down this safe column to the last row, which is guaranteed to succeed.
successful_attempts = 1

# The minimum value of n that guarantees a win is the total number of attempts
# required in this worst-case scenario.
n = max_failed_attempts + successful_attempts

# Print the explanation and the final equation for n.
print("To find the guaranteed number of attempts (n), we analyze the worst-case scenario.")
print(f"Number of columns on the board: {num_cols}")
print(f"Number of columns containing a monster (unsafe columns): {num_unsafe_cols}")
print("In the worst case, Turbo must test every unsafe column before finding the safe one.")
print(f"Number of attempts to identify all unsafe columns = {max_failed_attempts}")
print(f"Number of attempts for the final, successful path = {successful_attempts}")
print("\nThe final equation for n is the sum of the failed attempts and the one successful attempt:")
print(f"{max_failed_attempts} + {successful_attempts} = {n}")