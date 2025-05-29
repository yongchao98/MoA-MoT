def is_valid_combination(numbers, letters):
    # Constraints based on the feedback
    # Numbers must be 4 and 5
    if set(numbers) != {"4", "5"}:
        return False
    
    # Letters must be A and G in some order
    if set(letters) != {"A", "G"}:
        return False
    
    # Check specific constraints from each guess
    # Guess 4: 67AG - both numbers are incorrect and too large; both letters are correct but in the wrong positions
    if numbers == ["6", "7"] or numbers == ["7", "6"]:
        return False
    if letters == ["A", "G"]:
        return False
    
    # Guess 7: 35DU - one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect
    if numbers == ["3", "5"] or numbers == ["5", "3"]:
        return False
    if "D" in letters or "U" in letters:
        return False
    
    return True

def find_password():
    # Possible numbers and letters
    possible_numbers = ["4", "5"]
    possible_letters = ["A", "G"]
    
    # Try all combinations of numbers and letters
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for let1 in possible_letters:
                for let2 in possible_letters:
                    if let1 == let2:
                        continue
                    if is_valid_combination([num1, num2], [let1, let2]):
                        return [num1, num2, let1, let2]

# Find and print the password
password = find_password()
print(f"<<< {password} >>>")