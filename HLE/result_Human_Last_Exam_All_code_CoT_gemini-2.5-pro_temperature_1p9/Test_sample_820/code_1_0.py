import sys

# Increase recursion limit for deep calculations
sys.setrecursionlimit(2000)

# Initial configuration of disks on pegs
# Peg 0: [7, 3, 2], Peg 1: [1], Peg 2: [8, 6], Peg 3: [9, 5, 4], Peg 4: []
initial_pegs = {
    1: 1, 2: 0, 3: 0, 4: 3, 5: 3, 6: 2, 7: 0, 8: 2, 9: 3
}

# Memoization cache for our main solver
memo = {}
# List to store the explanation strings for final output
explanation_steps = []

def calculate_standard_moves():
    """Calculates the minimum moves for standard tower transfers."""
    T3 = [0] * 10
    T4 = [0] * 10
    T5 = [0] * 10

    for n in range(1, 10):
        # 3 pegs (standard Hanoi)
        T3[n] = 2 * T3[n - 1] + 1
        
        # 4 pegs (Frame-Stewart)
        T4[n] = min(2 * T4[n - j] + T3[j] for j in range(1, n))
        
        # 5 pegs (Frame-Stewart)
        T5[n] = min(2 * T5[n - k] + T4[k] for k in range(1, n))
        
    return T5

# Pre-calculate moves for a standard 5-peg Hanoi tower
T5 = calculate_standard_moves()

def solve(k, target_peg):
    """
    Recursively calculates the minimum moves to move disks 1..k to the target_peg.
    Returns a tuple (value, optimal_aux_peg_if_any).
    """
    if k == 0:
        return 0, None
        
    if (k, target_peg) in memo:
        return memo[(k, target_peg)]

    source_peg = initial_pegs[k]

    if source_peg == target_peg:
        # Disk k is already on the target peg for this subproblem.
        # Solve for the k-1 disks to be placed on top of it.
        res, _ = solve(k - 1, target_peg)
        optimal_aux = None
    else:
        # Disk k needs to be moved from source_peg to target_peg.
        # Find the best auxiliary peg to temporarily move disks 1..k-1.
        aux_pegs = [p for p in range(5) if p != source_peg and p != target_peg]
        
        min_sub_moves = float('inf')
        optimal_aux = -1

        for aux in aux_pegs:
            sub_moves, _ = solve(k - 1, aux)
            if sub_moves < min_sub_moves:
                min_sub_moves = sub_moves
                optimal_aux = aux
        
        # Total moves = move 1..k-1 to aux + move k to target + move 1..k-1 to target
        res = min_sub_moves + 1 + T5[k - 1]

    memo[(k, target_peg)] = (res, optimal_aux)
    return res, optimal_aux

def format_explanation(k, target_peg):
    """Generates the explanation string for a given step."""
    val, aux = memo.get((k, target_peg), (None, None))
    if val is None:
        return

    source_peg = initial_pegs[k]
    step_str = f"M({k}, {target_peg})"
    
    if source_peg == target_peg:
        # Disk k is already at the target for this subproblem
        prev_val, _ = memo.get((k - 1, target_peg), (0, None))
        explanation_steps.append(f"{step_str} = M({k-1}, {target_peg}) = {prev_val}")
    else:
        # Disk k needs to be moved
        sub_val, _ = memo.get((k - 1, aux), (0, None))
        tower_moves = T5[k-1]
        explanation_steps.append(f"{step_str} = M({k-1}, {aux}) + 1 + T5({k-1}) = {sub_val} + 1 + {tower_moves} = {val}")

if __name__ == "__main__":
    # Start the solver for all 9 disks, targeting Peg 4.
    total_moves, _ = solve(9, 4)

    # Generate explanations for each step by backtracking the solution
    for k in range(1, 10):
        for p in range(5):
            # We only need to format steps that were actually computed
            if (k, p) in memo:
                format_explanation(k, p)

    # Print the equation steps in a readable order
    print("Calculation Steps:")
    unique_steps = sorted(list(set(explanation_steps)), key=lambda x: (int(x.split('(')[1].split(',')[0]), x))
    for step in unique_steps:
        print(step)
    
    # Print the final answer
    print("\nFinal Answer:")
    print(f"The minimal amount of moves is {total_moves}.")
    print(f"<<<{total_moves}>>>")
