def check_numbers(n1, n2, guess_n1, guess_n2, feedback):
    correct_pos = 0
    wrong_pos = 0
    too_small = 0
    too_large = 0
    
    if n1 == guess_n1:
        correct_pos += 1
    elif n1 == guess_n2:
        wrong_pos += 1
    elif n1 < guess_n1:
        too_large += 1
    else:
        too_small += 1
        
    if n2 == guess_n2:
        correct_pos += 1
    elif n2 == guess_n1:
        wrong_pos += 1
    elif n2 < guess_n2:
        too_large += 1
    else:
        too_small += 1
        
    return correct_pos, wrong_pos, too_small, too_large

def check_letters(l1, l2, guess_l1, guess_l2, feedback):
    correct_pos = 0
    wrong_pos = 0
    too_early = 0
    too_late = 0
    
    if l1 == guess_l1:
        correct_pos += 1
    elif l1 == guess_l2:
        wrong_pos += 1
    elif l1 < guess_l1:
        too_late += 1
    else:
        too_early += 1
        
    if l2 == guess_l2:
        correct_pos += 1
    elif l2 == guess_l1:
        wrong_pos += 1
    elif l2 < guess_l2:
        too_late += 1
    else:
        too_early += 1
        
    return correct_pos, wrong_pos, too_early, too_late

def check_guess(password, guess, feedback):
    n1, n2 = int(password[0]), int(password[1])
    l1, l2 = password[2], password[3]
    
    guess_n1, guess_n2 = int(guess[0]), int(guess[1])
    guess_l1, guess_l2 = guess[2], guess[3]
    
    num_correct_pos, num_wrong_pos, num_too_small, num_too_large = check_numbers(n1, n2, guess_n1, guess_n2, feedback)
    let_correct_pos, let_wrong_pos, let_too_early, let_too_late = check_letters(l1, l2, guess_l1, guess_l2, feedback)
    
    if feedback == 1:
        return (num_wrong_pos == 1 and num_too_small == 1 and let_too_late == 2)
    elif feedback == 2:
        return (num_too_small == 2 and let_wrong_pos == 0 and let_correct_pos == 0)
    elif feedback == 3:
        return (num_correct_pos == 1 and num_too_large == 1 and 
                let_correct_pos == 1 and let_too_early == 1)
    else:  # feedback == 4
        return (num_correct_pos == 1 and num_too_small == 1 and 
                let_correct_pos == 1 and let_too_late == 1)

def find_password():
    numbers = list(range(10))
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    password = f"{n1}{n2}{l1}{l2}"
                    if (check_guess(password, "25YN", 1) and
                        check_guess(password, "42CQ", 2) and
                        check_guess(password, "76KC", 3) and
                        check_guess(password, "26SG", 4)):
                        return [str(n1), str(n2), l1, l2]

result = find_password()
print(result)