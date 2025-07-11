import fractions

# Part 1: How many pearls were there altogether?

# First, calculate the number of pearls that remained on the string.
# "seven shy of eleven times eleven"
remained_on_string = (11 * 11) - 7

# Let 'x' be the total number of pearls. The equation is:
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remained_on_string
# To solve for x, we rearrange the equation:
# x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = remained_on_string
# x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = remained_on_string
# We use the fractions module for precision.
sum_of_fractions = fractions.Fraction(1, 6) + fractions.Fraction(1, 5) + fractions.Fraction(1, 3) + fractions.Fraction(1, 10)
remaining_fraction = 1 - sum_of_fractions

# x = remained_on_string / remaining_fraction
total_pearls = remained_on_string / remaining_fraction

# Calculate the number of pearls for each part
on_floor = total_pearls * fractions.Fraction(1, 6)
on_bed = total_pearls * fractions.Fraction(1, 5)
saved_by_woman = total_pearls * fractions.Fraction(1, 3)
caught_by_lover = total_pearls * fractions.Fraction(1, 10)

print("--- Solving the riddle ---")
print(f"First, we solve for the total number of pearls.")
print(f"The equation representing the pearls is:")
print(f"Total Pearls = (Fell to floor) + (Fell on bed) + (Saved by woman) + (Caught by lover) + (Remained on string)")
print(f"{int(total_pearls)} = {int(on_floor)} + {int(on_bed)} + {int(saved_by_woman)} + {int(caught_by_lover)} + {int(remained_on_string)}\n")
print(f"There were {int(total_pearls)} pearls altogether.")

print("\n--- Answering the second question ---")
# Part 2: How many more are they gonna need?
fallen_pearls = total_pearls - remained_on_string
print(f"Total pearls that fell: {int(fallen_pearls)}")

found_pearls = fallen_pearls / 3
print(f"They managed to find back 1/3rd of them: {int(found_pearls)}")

needed_pearls = fallen_pearls - found_pearls
print(f"Therefore, the number of pearls they still need is: {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}")
print(f"\nThey are going to need {int(needed_pearls)} more pearls.")

# The final answer for the second question is provided below as requested.
# The format is <<<answer content>>>
# The answer is the number of pearls they still need.
final_answer = int(needed_pearls)