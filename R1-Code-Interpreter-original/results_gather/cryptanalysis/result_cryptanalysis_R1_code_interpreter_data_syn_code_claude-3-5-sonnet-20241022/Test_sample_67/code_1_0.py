def verify_guess(guess, target):
    # Returns feedback based on the guess compared to target
    feedback = []
    for i in range(4):
        if guess[i] == target[i]:
            feedback.append("correct position")
        elif guess[i] in target:
            feedback.append("wrong position")
        else:
            if i < 2:  # numbers
                if guess[i] > target[i]:
                    feedback.append("too large")
                else:
                    feedback.append("too small")
            else:  # letters
                if ord(guess[i]) > ord(target[i]):
                    feedback.append("too late")
                else:
                    feedback.append("too early")
    return feedback

def test_combination(combination):
    # Test against all known conditions
    test_cases = [
        ([7,9,'F','V'], ["too large", "too large", "incorrect", "incorrect"]),
        ([3,2,'P','Z'], ["too small", "too small", "too late", "too late"]),
        ([0,9,'E','F'], ["incorrect", "incorrect", "too early", "too early"]),
        ([5,8,'Q','D'], ["wrong position", "too large", "incorrect", "incorrect"]),
        ([7,9,'O','J'], ["too large", "too large", "incorrect", "incorrect"]),
        ([6,4,'L','Y'], ["too large", "wrong position", "correct position", "too late"]),
        ([4,8,'H','I'], ["correct position", "too large", "too early", "too early"]),
        ([4,5,'T','G'], ["correct position", "correct position", "incorrect", "incorrect"]),
        ([3,1,'I','B'], ["too small", "too small", "too early", "too early"]),
        ([9,4,'V','W'], ["too large", "wrong position", "too late", "too late"]),
        ([7,0,'X','N'], ["incorrect", "incorrect", "too late", "too late"]),
        ([7,0,'B','I'], ["incorrect", "incorrect", "too early", "too early"]),
        ([8,9,'U','G'], ["too large", "too large", "incorrect", "incorrect"]),
        ([7,0,'K','G'], ["incorrect", "incorrect", "too early", "too early"])
    ]
    
    for test_case, expected_feedback in test_cases:
        feedback = verify_guess(test_case, combination)
        if feedback != expected_feedback:
            return False
    return True

# Test the proposed solution [4,5,'L','M']
solution = [4,5,'L','M']
print(f"Testing solution {solution}")
print(f"Solution is valid: {test_combination(solution)}")

# Let's also verify there are no other valid solutions
other_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:
                    continue
                candidate = [n1, n2, l1, l2]
                if candidate != solution and test_combination(candidate):
                    other_solutions.append(candidate)

print(f"Number of other valid solutions found: {len(other_solutions)}")
if other_solutions:
    print("Other solutions:", other_solutions)