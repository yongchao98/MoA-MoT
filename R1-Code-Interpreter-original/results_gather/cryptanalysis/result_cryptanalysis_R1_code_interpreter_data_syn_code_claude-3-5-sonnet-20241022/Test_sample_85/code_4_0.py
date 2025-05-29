def verify_guess(guess):
    n1, n2, l1, l2 = guess
    
    # Helper function to check if numbers match feedback
    def check_numbers(test_nums, feedback):
        g_n1, g_n2 = test_nums
        if feedback == "both_incorrect":
            return n1 != g_n1 and n2 != g_n2
        elif feedback == "both_incorrect_large":
            return n1 != g_n1 and n2 != g_n2 and n1 < g_n1 and n2 < g_n2
        elif feedback == "one_correct_position":
            return (n1 == g_n1 and n2 != g_n2) or (n2 == g_n2 and n1 != g_n1)
        elif feedback == "one_wrong_position":
            return (n1 == g_n2 or n2 == g_n1) and n1 != g_n1 and n2 != g_n2
    
    # Helper function to check if letters match feedback
    def check_letters(test_letters, feedback):
        g_l1, g_l2 = test_letters
        if feedback == "both_incorrect":
            return l1 != g_l1 and l2 != g_l2
        elif feedback == "one_correct_position":
            return (l1 == g_l1 and l2 != g_l2) or (l2 == g_l2 and l1 != g_l1)
        elif feedback == "one_wrong_position":
            return (l1 == g_l2 or l2 == g_l1) and l1 != g_l1 and l2 != g_l2

    # Check all conditions
    conditions = [
        ((8,6,'A','F'), ("both_incorrect", "both_incorrect")),
        ((9,8,'L','F'), ("both_incorrect_large", "both_incorrect")),
        ((2,0,'X','K'), ("one_correct_position", "both_incorrect")),
        ((3,8,'T','M'), ("both_incorrect", "both_incorrect")),
        ((8,1,'J','K'), ("both_incorrect", "both_incorrect")),
        ((2,3,'V','Q'), ("one_correct_position", "both_incorrect")),
        ((8,5,'I','C'), ("both_incorrect", "both_incorrect")),
        ((0,3,'O','E'), ("both_incorrect", "both_incorrect")),
        ((7,8,'W','K'), ("one_wrong_position", "one_correct_position")),
        ((6,8,'K','C'), ("both_incorrect", "both_incorrect")),
        ((4,7,'O','M'), ("one_correct_position", "both_incorrect")),
        ((5,8,'Z','Q'), ("both_incorrect", "both_incorrect")),
        ((4,9,'D','W'), ("both_incorrect", "one_wrong_position"))
    ]
    
    for (test_n1, test_n2, test_l1, test_l2), (num_feedback, let_feedback) in conditions:
        if not check_numbers((test_n1, test_n2), num_feedback):
            return False
        if not check_letters((test_l1, test_l2), let_feedback):
            return False
    
    return True

# Test all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:  # Numbers can't repeat
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:  # Letters can't repeat
                    continue
                guess = (n1, n2, l1, l2)
                if verify_guess(guess):
                    valid_solutions.append(guess)

print(f"Found {len(valid_solutions)} solution(s):")
for solution in valid_solutions:
    print(solution)