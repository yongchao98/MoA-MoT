import random

def illustrate_non_positive_recurrence():
    """
    Illustrates that a chain satisfying the given conditions is not
    positive recurrent using the example of a Symmetric Random Walk (SRW) on Z.

    The SRW is defined by p(x, x+1) = p(x, x-1) = 0.5. It is irreducible.
    Let A = {0}, which is a finite set.
    Let f(x) = x^2. This function is non-negative and f(x) -> infinity as |x| -> infinity.

    We first verify the condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A.
    Then, we simulate the chain to show that the mean return time to the origin
    is not finite, a hallmark of non-positive recurrent chains.
    """

    print("--- Part 1: Verifying the condition for the SRW example ---")
    
    # Let's verify the condition for a specific x, e.g., x=5.
    x = 5
    p = 0.5
    fx = x**2
    fx_plus_1 = (x + 1)**2
    fx_minus_1 = (x - 1)**2
    expected_f = p * fx_plus_1 + p * fx_minus_1
    drift = expected_f - fx

    print(f"Let's check the condition for x = {x}, with f(x) = x^2 and A = {{0}}.")
    print(f"The condition is: sum_y(p(x,y)f(y)) - f(x) >= 0")
    print(f"For the SRW, p({x},{x+1}) = {p} and p({x},{x-1}) = {p}.")
    print("\nCalculating each term of the equation:")
    print(f"f({x}) = {x}^2 = {fx}")
    print(f"E[f(X_1)|X_0={x}] = {p}*f({x+1}) + {p}*f({x-1})")
    print(f"               = {p}*({x+1})^2 + {p}*({x-1})^2")
    print(f"               = {p}*{fx_plus_1} + {p}*{fx_minus_1}")
    print(f"               = {p*fx_plus_1} + {p*fx_minus_1} = {expected_f}")
    print(f"\nThe final result of the equation is:")
    print(f"E[f(X_1)|X_0={x}] - f(x) = {expected_f} - {fx} = {drift}")
    print(f"Since {drift} >= 0, the condition is met for x={x}. It holds for all x != 0.")

    print("\n--- Part 2: Simulating the SRW to estimate mean return time ---")
    
    num_simulations = 10000
    max_steps_per_simulation = 50000 # A cap to prevent infinite loops
    
    return_times = []
    failed_to_return = 0

    for _ in range(num_simulations):
        position = 0
        # Start one step away from the origin to measure return time
        position = 1 if random.random() < 0.5 else -1
        
        for steps in range(1, max_steps_per_simulation + 1):
            # Move left or right with equal probability
            position += 1 if random.random() < 0.5 else -1
            
            if position == 0:
                return_times.append(steps)
                break
        else: # This else belongs to the for loop, executed if the loop finishes without break
            failed_to_return += 1

    if return_times:
        average_return_time = sum(return_times) / len(return_times)
        print(f"Number of simulations: {num_simulations}")
        print(f"Max steps per simulation: {max_steps_per_simulation}")
        print(f"Number of simulations where the chain returned to 0: {len(return_times)}")
        print(f"Number of simulations where the chain did not return: {failed_to_return}")
        print(f"Average return time (for successful returns): {average_return_time:.2f}")
    else:
        print("The chain did not return to the origin in any simulation.")

    print("\nInterpretation:")
    print("For a positive recurrent chain, the mean return time is finite. Our simulation shows a very large average return time,")
    print("and a significant number of paths failed to return within the step limit. This suggests the mean return time is infinite,")
    print("which is characteristic of a null recurrent or transient chain. Thus, the chain is not positive recurrent.")

if __name__ == '__main__':
    illustrate_non_positive_recurrence()