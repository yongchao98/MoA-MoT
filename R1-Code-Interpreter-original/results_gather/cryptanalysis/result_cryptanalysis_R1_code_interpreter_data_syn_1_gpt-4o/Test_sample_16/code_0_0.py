# Possible numbers and letters based on deductions
possible_numbers = [8, 9]
possible_letters = ['N', 'A']

# Check each condition
def check_conditions(numbers, letters):
    # Condition 1: 80KR
    if (numbers[0] == 8 or numbers[1] == 8) and (numbers[0] != 0 and numbers[1] != 0) and ('K' not in letters and 'R' not in letters):
        # Condition 2: 87NP
        if (numbers[0] == 8 or numbers[1] == 8) and ('N' in letters) and ('P' not in letters):
            # Condition 3: 47AZ
            if ('A' in letters) and ('Z' not in letters):
                # Condition 4: 15RF
                if (numbers[0] > 5 and numbers[1] > 5) and ('R' not in letters and 'F' not in letters):
                    return True
    return False

# Check the combination
if check_conditions(possible_numbers, possible_letters):
    print(["8", "9", "N", "A"])
else:
    print("No valid combination found.")