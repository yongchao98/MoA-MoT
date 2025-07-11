import numpy as np
from scipy.optimize import linprog

def check_learnability(operator_name, truth_table):
    """
    Checks if a 2-variable boolean operator is learnable by the additively
    separable model. It does this by checking if a system of linear
    inequalities has a feasible solution.

    The model's score function S(p,q) must satisfy:
    S(p,q) = F(p) + G(q). Let F(p)=x_p*p, G(q)=x_q*q and the bias be b.
    S(0,0) = b
    S(0,1) = x_q + b
    S(1,0) = x_p + b
    S(1,1) = x_p + x_q + b

    We are looking for a solution for variables [x_p, x_q, b].
    """
    # The linear system is defined by A_ub * x <= b_ub
    # We formulate our problem as A * v >= epsilon, where v = [x_p, x_q, b].
    # This is equivalent to -A * v <= -epsilon.
    A = []
    # truth_table is for inputs (0,0), (0,1), (1,0), (1,1)
    targets = {
        (0, 0): truth_table[0],
        (0, 1): truth_table[1],
        (1, 0): truth_table[2],
        (1, 1): truth_table[3],
    }
    
    # Coefficients for [x_p, x_q, b]
    scores_coeffs = {
        (0, 0): [0, 0, 1],
        (0, 1): [0, 1, 1],
        (1, 0): [1, 0, 1],
        (1, 1): [1, 1, 1],
    }

    for inputs, target in targets.items():
        coeffs = np.array(scores_coeffs[inputs])
        if target == 1:
            # We need Score > 0, which means -Score < 0.
            A.append(-coeffs)
        else: # target == 0
            # We need Score < 0.
            A.append(coeffs)
            
    A_ub = np.array(A)
    # We want A*v > 0, which we can write as A*v >= 1 for simplicity
    # (any separating plane can be scaled).
    # So, -A*v <= -1.
    b_ub = -np.ones(4)

    # We want to find if a feasible solution exists. The objective function doesn't matter.
    c = np.zeros(3)

    # Solve the linear program
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(None, None))

    return res.success

def main():
    """
    Main function to define operators and test their learnability.
    """
    # Define the truth tables for mixing-dimension operators
    # The order is for (p1,q2) inputs: (0,0), (0,1), (1,0), (1,1)
    operators = {
        "X' (XOR)":          [0, 1, 1, 0],
        "C' (Conjunction)":  [0, 0, 0, 1],
        "D' (Disjunction)":  [0, 1, 1, 1],
        "E' (Equivalence)":  [1, 0, 0, 1],
        "I' (Implication)":  [1, 1, 0, 1],
    }

    print("Analyzing learnability of mixing-dimension operators...")
    
    unlearnable_operators = []
    for name, truth_table in operators.items():
        is_learnable = check_learnability(name, truth_table)
        print(f"- Operator {name}: {'Learnable' if is_learnable else 'NOT Learnable'}")
        if not is_learnable:
            unlearnable_operators.append(name.split(" ")[0])
    
    print("\n--- Conclusion ---")
    print("The operators that can not be learned with the heuristic representation are:")
    print(' '.join(unlearnable_operators))

if __name__ == "__main__":
    main()