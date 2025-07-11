from fractions import Fraction

# --- Part 1: Find the total number of pearls ---

# Calculate the number of pearls remaining on the string
# "a seven shy of eleven times eleven"
remaining_on_string = (11 * 11) - 7

# Define the fractions of the total pearls that fell
fractions_fallen = [Fraction(1, 6), Fraction(1, 5), Fraction(1, 3), Fraction(1, 10)]
sum_of_fallen_fractions = sum(fractions_fallen)

# The fraction of pearls remaining on the string is 1 minus the sum of the fallen fractions
fraction_on_string = 1 - sum_of_fallen_fractions

# The total number of pearls (x) can be found with the equation:
# x * fraction_on_string = remaining_on_string
# So, x = remaining_on_string / fraction_on_string
total_pearls = int(remaining_on_string / fraction_on_string)

# --- Print the results for Part 1 ---
print("To find the total number of pearls, we solve the following equation:")

# Calculate the number of pearls for each fraction for the printout
parts = [int(total_pearls * f) for f in fractions_fallen]
print(f"{total_pearls} = {parts[0]} (fell to floor) + {parts[1]} (fell on bed) + {parts[2]} (woman saved) + {parts[3]} (lover caught) + {remaining_on_string} (remained on string)")
print("-" * 20)
print(f"There were {total_pearls} pearls altogether.")
print("-" * 20)

# --- Part 2: Find how many more pearls they need ---

# Calculate the total number of fallen pearls
fallen_pearls = total_pearls - remaining_on_string

# They find 1/3 of the fallen pearls, so they still need the other 2/3
needed_pearls = int(fallen_pearls * Fraction(2, 3))

# --- Print the results for Part 2 ---
print(f"Total fallen pearls: {fallen_pearls}")
print(f"If they find 1/3 of the fallen pearls, they are still missing the other 2/3.")
print(f"Therefore, they are going to need {needed_pearls} more pearls.")

# This is a hidden value for the final answer format
final_answer = f"{total_pearls}, {needed_pearls}"
# print(f"<<<{final_answer}>>>")