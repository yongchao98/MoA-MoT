def check_guess(numbers, letters, guess_num, guess_let, feedback):
    # Check if a guess matches the feedback
    if len(set(numbers)) != 2 or len(set(letters)) != 2:  # Check for repeats
        return False
    
    # Convert to lists for position comparison
    num_guess = list(str(guess_num).zfill(2))
    let_guess = list(guess_let)
    num_actual = list(map(str, numbers))
    
    # Check number conditions
    if feedback == "both_small":
        return all(int(n) > int(g) for n, g in zip(num_actual, num_guess))
    elif feedback == "one_correct_wrong_pos_one_small":
        correct_but_wrong_pos = any(str(n) in num_guess and str(n) != num_guess[i] 
                                  for i, n in enumerate(numbers))
        one_too_small = any(n > int(g) for n, g in zip(numbers, num_guess))
        return correct_but_wrong_pos and one_too_small
    elif feedback == "both_incorrect":
        return all(str(n) not in num_guess for n in numbers)
    
    return True

def check_letters(numbers, letters, guess_letters, feedback):
    if feedback == "one_correct_wrong_pos_one_early":
        return (any(l in guess_letters and l != guess_letters[i] for i, l in enumerate(letters)) and
                any(ord(l) > ord(g) for l, g in zip(letters, guess_letters)))
    elif feedback == "both_early":
        return all(ord(l) > ord(g) for l, g in zip(letters, guess_letters))
    elif feedback == "both_incorrect":
        return all(l not in guess_letters for l in letters)
    return True

def find_password():
    # Try all possible combinations
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:
                continue
            numbers = (n1, n2)
            
            # Check number conditions
            if not all([
                check_guess(numbers, None, 30, None, "both_small"),
                check_guess(numbers, None, 60, None, "one_correct_wrong_pos_one_small"),
                check_guess(numbers, None, 23, None, "both_small"),
                check_guess(numbers, None, 03, None, "both_small"),
                check_guess(numbers, None, 79, None, "both_incorrect")
            ]):
                continue
                
            # If numbers pass, try letters
            for l1 in range(ord('A'), ord('Z')+1):
                for l2 in range(ord('A'), ord('Z')+1):
                    if l1 == l2:
                        continue
                    letters = (chr(l1), chr(l2))
                    
                    # Check letter conditions
                    if all([
                        check_letters(numbers, letters, "AU", "one_correct_wrong_pos_one_early"),
                        check_letters(numbers, letters, "KI", "both_early"),
                        check_letters(numbers, letters, "QB", "both_incorrect"),
                        check_letters(numbers, letters, "HC", "both_early"),
                        check_letters(numbers, letters, "MI", "both_incorrect")
                    ]):
                        return [str(n1), str(n2), letters[0], letters[1]]
    
    return None

result = find_password()
print(result)