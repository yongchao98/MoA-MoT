import random
import math

def calculate_conflicts(board, n):
    conflicts = 0
    for i in range(n):
        for j in range(i + 1, n):
            if board[i] == board[j] or abs(board[i] - board[j]) == abs(i - j):
                conflicts += 1
    return conflicts

def simulated_annealing(n, x_constraints):
    # Initial random placement of queens
    board = [-1] * n
    for i in range(n):
        if i not in x_constraints:
            board[i] = random.choice([j for j in range(n) if j not in x_constraints[i]])

    temperature = 4000
    cooling_rate = 0.99
    current_conflicts = calculate_conflicts(board, n)

    while temperature > 0.1 and current_conflicts > 0:
        new_board = board[:]
        row = random.choice([i for i in range(n) if i not in x_constraints])
        new_col = random.choice([j for j in range(n) if j not in x_constraints[row]])
        new_board[row] = new_col

        new_conflicts = calculate_conflicts(new_board, n)
        delta = new_conflicts - current_conflicts

        if delta < 0 or random.uniform(0, 1) < math.exp(-delta / temperature):
            board = new_board
            current_conflicts = new_conflicts

        temperature *= cooling_rate

    if current_conflicts == 0:
        return board
    else:
        return None

def place_queens_with_constraints():
    n = 8
    x_constraints = {
        2: {2},
        3: {1, 4},
        7: {6}
    }

    solution = simulated_annealing(n, x_constraints)
    if solution:
        print("<<<" + ", ".join(f"{i} {solution[i]}" for i in range(n)) + ">>>")
    else:
        print("No solution found")

place_queens_with_constraints()