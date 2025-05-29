from itertools import permutations

# Define possible numbers and letters
numbers = set(range(10))
letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Define constraints based on feedback
constraints = [
    (set(), set("SW"), "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    (set(), set("DW"), "both numbers are incorrect; both letters are incorrect."),
    (set(), set("ID"), "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ({6}, set("MG"), "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect."),
    (set(), set("RQ"), "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    ({3}, set("CM"), "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect."),
    (set(), set("YF"), "both numbers are incorrect; both letters are incorrect."),
    (set(), set("YU"), "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    ({3}, set("AG"), "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet."),
    ({4}, set("EB"), "one number is correct but in the wrong position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet."),
    ({3}, set("OB"), "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect."),
    ({3}, set("LK"), "one number is correct and in the correct position; one number is incorrect and too large; one letter is correct but in the wrong position; one letter is incorrect and incorrect."),
    (set(), set("UW"), "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    (set(), set("EU"), "both numbers are incorrect; both letters are incorrect."),
    ({3}, set("GE"), "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ({6}, set("VO"), "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect and too late in the alphabet."),
    (set(), set("IK"), "both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet.")
]

# Filter possible numbers and letters based on constraints
possible_numbers = numbers - {0, 1, 2, 4, 5, 6, 7, 9}
possible_letters = letters - set("SWDWRQCMYFUE")

# Check all permutations of possible numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        # Check if the current combination satisfies all constraints
        valid = True
        for correct_nums, incorrect_letters, feedback in constraints:
            if "both numbers are incorrect" in feedback:
                if any(n in num_perm for n in correct_nums):
                    valid = False
                    break
            if "one number is correct and in the correct position" in feedback:
                if not any(num_perm[i] in correct_nums for i in range(2)):
                    valid = False
                    break
            if "one number is correct but in the wrong position" in feedback:
                if not any(num_perm[i] in correct_nums for i in range(2)):
                    valid = False
                    break
            if "both letters are incorrect" in feedback:
                if any(l in letter_perm for l in incorrect_letters):
                    valid = False
                    break
            if "one letter is correct but in the wrong position" in feedback:
                if not any(letter_perm[i] in incorrect_letters for i in range(2)):
                    valid = False
                    break
        if valid:
            # Output the valid combination
            print(f"<<< {list(num_perm) + list(letter_perm)} >>>")
            break