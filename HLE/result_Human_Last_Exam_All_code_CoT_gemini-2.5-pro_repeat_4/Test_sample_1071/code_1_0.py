import sympy

def solve_for_t_bounds():
    """
    This function programmatically solves for the bounds of t based on the problem's derivation.
    """
    print("This script finds the lower and upper bounds of t based on the given requirements.")
    
    # 1. Define the variable and initial constraints
    t = sympy.Symbol('t', real=True)
    
    # Constraint from problem statement: interval is [-1, t], so t >= -1
    constraint1 = sympy.Ge(t, -1)
    
    # Constraint to avoid division by zero: a0+a2 != 0.
    # The sum range is [-2, 2*t]. For 0 not to be in this range, we need 2*t < 0.
    constraint2 = sympy.Lt(t, 0)
    print(f"\nStep 1: To prevent (a0 + a2) from being 0, we must have t < 0.")

    # 2. Define the main condition based on set inclusion.
    # The range of sums S is [-2, 2*t].
    # The set of reciprocals {1/s | s in S} is [1/(2*t), -1/2].
    # This set must be a subset of S. This leads to two inequalities.
    
    # Inequality 1: -1/2 <= 2*t
    inequality1 = sympy.Ge(2*t, -sympy.Rational(1, 2))
    
    # Inequality 2: -2 <= 1/(2*t)
    # Since t < 0, multiplying by 2*t reverses the inequality: -4*t >= 1
    inequality2 = sympy.Ge(-4*t, 1)

    print(f"Step 2: The condition implies two inequalities for t:")
    print(f"  a) -1/2 <= 2*t  (derived from the upper bounds)")
    print(f"  b) -2 <= 1/(2*t) (derived from the lower bounds)")

    # 3. Solve the system of inequalities
    solution_set = sympy.solveset(
        sympy.And(constraint1, constraint2, inequality1, inequality2), 
        t, 
        domain=sympy.S.Reals
    )
    
    print("\nStep 3: Solving these inequalities simultaneously yields a single value for t.")

    if isinstance(solution_set, sympy.FiniteSet) and len(solution_set) == 1:
        solution_value = list(solution_set)[0]
        lower_bound = float(solution_value)
        upper_bound = float(solution_value)
        
        print(f"\nThe only value of t that satisfies all conditions is t = {solution_value}.")
        print("Therefore, the lower and upper bounds are the same.")
        
        print("\nLower and Upper Bounds:")
        # The prompt asks to output the numbers for the final answer
        print(f"{lower_bound} {upper_bound}")
        return lower_bound, upper_bound
    else:
        print("Could not find a unique solution. The solution set is:", solution_set)
        return None, None

solve_for_t_bounds()