import math

def analyze_strategies():
    """
    Analyzes the optimal resource distribution strategy for a mother bird
    by testing various survival functions against Fair and Unfair strategies.
    """

    # 1. Define problem parameters
    n = 5       # Number of offspring
    R = 10.0    # Total resources
    r_max = 4.0   # Max resources per offspring

    print(f"Problem Setup: n={n}, R={R}, r_max={r_max}\n")

    # 2. Implement the strategies
    def get_fair_allocation(n, R, r_max):
        """Calculates the allocation for the Fair Strategy."""
        # Check if feasible
        if R / n > r_max:
            return None
        return [R / n] * n

    def get_unfair_allocation(n, R, r_max):
        """Calculates the allocation for the Unfair Strategy."""
        allocation = [0.0] * n
        resources_left = R
        k = int(resources_left / r_max)
        
        # Give r_max to k offspring
        for i in range(k):
            if resources_left >= r_max:
                allocation[i] = r_max
                resources_left -= r_max
        
        # Give remainder to the next offspring
        if k < n and resources_left > 0:
            allocation[k] = resources_left
        
        return allocation

    fair_alloc = get_fair_allocation(n, R, r_max)
    unfair_alloc = get_unfair_allocation(n, R, r_max)

    print(f"Fair Strategy Allocation: {fair_alloc}")
    print(f"Unfair Strategy Allocation: {unfair_alloc}\n")
    print("---" * 15)

    def evaluate(s_func, s_name, alloc, alloc_name):
        """Evaluates a strategy and prints the detailed result."""
        results = [s_func(r) for r in alloc]
        total_survival = sum(results)
        
        # Format the equation string
        equation_parts = [f"s({r:.1f})" for r in alloc]
        equation_str = " + ".join(equation_parts)
        
        value_parts = [f"{res:.2f}" for res in results]
        value_str = " + ".join(value_parts)
        
        print(f"{alloc_name} for s(r)={s_name}:")
        print(f"  Sum = {equation_str}")
        print(f"      = {value_str}")
        print(f"      = {total_survival:.2f}\n")
        return total_survival

    # 3. Analyze each statement with example functions

    # Statement 1: "If s is strictly increasing, then the fair strategy is always optimal."
    print("Analysis for Statement 1:")
    # We test with a strictly increasing AND convex function.
    s1 = lambda r: r**2
    s1_name = "r^2 (Increasing, Convex)"
    
    fair_score = evaluate(s1, s1_name, fair_alloc, "Fair")
    unfair_score = evaluate(s1, s1_name, unfair_alloc, "Unfair")
    
    print(f"Result: Unfair ({unfair_score:.2f}) > Fair ({fair_score:.2f}).")
    print("Conclusion: The Unfair strategy is better here. Therefore, Statement 1 is FALSE.\n")
    print("---" * 15)


    # Statement 2: "If s is strictly decreasing, then the unfair strategy is always optimal."
    print("Analysis for Statement 2:")
    # We test with a strictly decreasing AND concave function.
    s2 = lambda r: -r**2
    s2_name = "-r^2 (Decreasing, Concave)"

    fair_score = evaluate(s2, s2_name, fair_alloc, "Fair")
    unfair_score = evaluate(s2, s2_name, unfair_alloc, "Unfair")

    print(f"Result: Fair ({fair_score:.2f}) > Unfair ({unfair_score:.2f}).")
    print("Conclusion: The Fair strategy is better here. Therefore, Statement 2 is FALSE.\n")
    print("---" * 15)

    
    # Statement 3: "If s is concave increasing, then fair is optimal. If s is concave decreasing, then unfair is optimal."
    print("Analysis for Statement 3:")
    print("Part A: Testing 'concave increasing' -> fair is optimal")
    # Using sqrt(r), which is concave and increasing
    s3 = lambda r: math.sqrt(r) if r >= 0 else 0
    s3_name = "sqrt(r) (Increasing, Concave)"
    fair_score = evaluate(s3, s3_name, fair_alloc, "Fair")
    unfair_score = evaluate(s3, s3_name, unfair_alloc, "Unfair")
    print(f"Result: Fair ({fair_score:.2f}) > Unfair ({unfair_score:.2f}). This part holds.")

    print("\nPart B: Testing 'concave decreasing' -> unfair is optimal")
    # Using the same function as for Statement 2 (-r^2), which is concave and decreasing
    fair_score = evaluate(s2, s2_name, fair_alloc, "Fair")
    unfair_score = evaluate(s2, s2_name, unfair_alloc, "Unfair")
    print(f"Result: Fair ({fair_score:.2f}) > Unfair ({unfair_score:.2f}). The Fair strategy is optimal.")
    print("Conclusion: Part B of the statement is false. Therefore, Statement 3 is FALSE.\n")
    print("---" * 15)


    # Statement 4: "If s is concave then the fair strategy is always optimal, regardless of whether it is increasing or decreasing."
    print("Analysis for Statement 4:")
    print("This statement claims fair is always optimal for any concave function.")
    print("Our test for s(r)=sqrt(r) (concave, increasing) showed fair was optimal.")
    print("Our test for s(r)=-r^2 (concave, decreasing) showed fair was optimal.")
    print("This aligns with Jensen's inequality, which proves that for a concave function, the average of values is greater than or equal to the value of the average. So n * s(R/n) >= sum(s(r_i)).")
    print("Conclusion: Statement 4 is TRUE.\n")
    print("---" * 15)

    # Final Conclusion
    print("\nFINAL SUMMARY:")
    print("Statement 1: FALSE")
    print("Statement 2: FALSE")
    print("Statement 3: FALSE")
    print("Statement 4: TRUE")
    print("\nThe only correct statement is [4].")

if __name__ == '__main__':
    analyze_strategies()