def solve_jug_problem_efficient():
    moves = [
        'fill B',    # Step 1: Fill Jug B (11 liters)
        'pour B->A', # Step 2: Pour from B to A (A: 6, B: 5)
        'empty A',   # Step 3: Empty Jug A (A: 0, B: 5)
        'pour B->A', # Step 4: Pour from B to A (A: 5, B: 0)
        'fill B',    # Step 5: Fill Jug B (A: 5, B: 11)
        'pour B->A'  # Step 6: Pour from B to A (A: 6, B: 3)
    ]
    return moves

solution_efficient = solve_jug_problem_efficient()
print(solution_efficient)