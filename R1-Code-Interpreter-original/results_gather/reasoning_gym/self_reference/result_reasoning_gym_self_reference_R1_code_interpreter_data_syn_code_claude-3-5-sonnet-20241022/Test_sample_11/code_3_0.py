def analyze_solution(true_statements):
    statements = [False] * 7
    for i in true_statements:
        statements[i-1] = True
    
    num_true = len(true_statements)
    num_false = 7 - num_true
    
    print(f"\nAnalyzing solution with true statements: {true_statements}")
    print(f"Number of true statements: {num_true}")
    print(f"Number of false statements: {num_false}")
    
    # Analyze each statement
    print("\nStatement analysis:")
    
    # Statement 1: "At least 4 are true"
    print(f"1. Claims 'at least 4 true': {num_true >= 4} (Statement is {statements[0]})")
    
    # Statement 2: "At most 5 are false"
    print(f"2. Claims 'at most 5 false': {num_false <= 5} (Statement is {statements[1]})")
    
    # Statement 3: "Exactly 7 are true"
    print(f"3. Claims 'exactly 7 true': {num_true == 7} (Statement is {statements[2]})")
    
    # Statement 4: "Exactly 1 is false"
    print(f"4. Claims 'exactly 1 false': {num_false == 1} (Statement is {statements[3]})")
    
    # Statement 5: "Statement 3 XOR Statement 4"
    xor_result = (statements[2] != statements[3]) and (statements[2] or statements[3])
    print(f"5. Claims '3 XOR 4': {xor_result} (Statement is {statements[4]})")
    
    # Statement 6: "Number of true statements is prime"
    print(f"6. Claims 'true count is prime': {num_true in [2,3,5,7]} (Statement is {statements[5]})")
    
    # Statement 7: "Number of false statements is composite"
    print(f"7. Claims 'false count is composite': {num_false in [4,6]} (Statement is {statements[6]})")

# Analyze each solution
solutions = [[], [2,6], [7], [2,6,7]]
for sol in solutions:
    analyze_solution(sol)