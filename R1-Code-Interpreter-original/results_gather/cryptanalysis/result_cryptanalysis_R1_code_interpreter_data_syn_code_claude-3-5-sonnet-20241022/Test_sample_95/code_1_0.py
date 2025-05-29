def verify_solution(solution):
    def check_numbers(guess_nums, feedback):
        if feedback == "both numbers are incorrect":
            return all(n not in solution[:2] for n in guess_nums)
        elif feedback == "one number is correct but in the wrong position":
            return any(n in solution[:2] and guess_nums.index(n) != solution.index(n) for n in guess_nums)
        elif feedback == "one number is correct and in the correct position":
            return any(n == solution[i] for i, n in enumerate(guess_nums))
        return True

    def check_letters(guess_letters, feedback):
        if feedback == "both letters are incorrect":
            return all(l not in solution[2:] for l in guess_letters)
        elif feedback == "one letter is correct and in the correct position":
            return any(l == solution[i+2] for i, l in enumerate(guess_letters))
        return True

    conditions = [
        ("79", "AG", "one number is correct but in the wrong position", "both letters are incorrect"),
        ("15", "BN", "one number is correct and in the correct position", "both letters are incorrect"),
        ("46", "QN", "both numbers are incorrect", "both letters are incorrect"),
        ("16", "BC", "both numbers are incorrect", "both letters are incorrect"),
        ("81", "BA", "both numbers are incorrect", "both letters are incorrect"),
        ("69", "IH", "one number is correct but in the wrong position", "both letters are incorrect"),
        ("05", "WO", "one number is correct and in the correct position", "both letters are incorrect"),
        ("74", "PV", "both numbers are incorrect", "one letter is correct and in the correct position"),
        ("26", "BF", "both numbers are incorrect", "both letters are incorrect"),
        ("63", "TH", "both numbers are incorrect", "both letters are incorrect"),
        ("74", "HT", "both numbers are incorrect", "both letters are incorrect"),
        ("63", "OY", "both numbers are incorrect", "both letters are incorrect"),
        ("06", "XR", "both numbers are incorrect", "both letters are incorrect"),
        ("63", "HS", "both numbers are incorrect", "both letters are incorrect"),
        ("13", "EX", "both numbers are incorrect and too small", "one letter is correct and in the correct position")
    ]

    for guess_nums, guess_letters, num_feedback, letter_feedback in conditions:
        if not (check_numbers(guess_nums, num_feedback) and check_letters(guess_letters, letter_feedback)):
            print(f"Failed condition: {guess_nums}{guess_letters}")
            return False
    return True

solution = "95EV"
print(f"Solution {solution} is {'valid' if verify_solution(solution) else 'invalid'}")