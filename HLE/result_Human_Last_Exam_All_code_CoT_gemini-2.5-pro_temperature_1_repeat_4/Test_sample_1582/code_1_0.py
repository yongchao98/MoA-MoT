import random
import matplotlib.pyplot as plt

def demonstrate_non_positive_recurrence():
    """
    This script demonstrates that if a Markov chain has a function f with certain properties,
    it cannot be positive recurrent.

    We use the example of a Symmetric Random Walk on Integers (Z).
    - State Space Sigma = Z = {..., -2, -1, 0, 1, 2, ...}
    - Transition probabilities: p(x, x+1) = 0.5, p(x, x-1) = 0.5.
    - This chain is known to be null recurrent, which is a type of non-positive recurrence.

    We need to find a function f and a finite set A that satisfy the user's conditions.
    - Let f(x) = x^2. This function is non-negative and f(x) -> infinity as |x| -> infinity.
    - Let the finite set A be the empty set, A = {}.

    Now, we check the main condition: is the drift non-negative for all x not in A?
    Drift(x) = E[f(X_n+1) | X_n=x] - f(x) >= 0 for all x in Z.
    """
    print("--- Verifying the Conditions for the Example ---")
    print("Markov Chain: Symmetric Random Walk on Integers")
    print("Function f(x) = x^2")
    print("Finite Set A = {} (the empty set)")
    print("\nChecking the drift condition: E[f(X_n+1) | X_n=x] - f(x) >= 0")

    # Let's verify the drift calculation for an example point, x=10
    x = 10
    p_right = 0.5
    p_left = 0.5
    
    # The next state can be x+1 or x-1
    x_plus_1 = x + 1
    x_minus_1 = x - 1

    # f(x) = x^2
    f = lambda val: val**2
    
    # The drift calculation
    expected_f_next = p_right * f(x_plus_1) + p_left * f(x_minus_1)
    current_f = f(x)
    drift = expected_f_next - current_f

    print(f"\nFor x = {x}:")
    print(f"E[f(X_n+1)|X_n={x}] = p(x,x+1)*f(x+1) + p(x,x-1)*f(x-1)")
    # The prompt requested to output each number in the final equation
    print("Final Equation (Calculation):")
    print(f"{p_right} * {x_plus_1}^2 + {p_left} * {x_minus_1}^2 - {x}^2")
    print(f"= {p_right} * {f(x_plus_1)} + {p_left} * {f(x_minus_1)} - {f(x)}")
    print(f"= {p_right * f(x_plus_1)} + {p_left * f(x_minus_1)} - {f(x)}")
    print(f"= {expected_f_next} - {current_f} = {drift}")

    print(f"\nThe calculated drift is {drift}, which is >= 0. This holds for any x.")
    print("The conditions are satisfied. The theory implies the chain is not positive recurrent.")
    print("\n--- Simulating the Chain to Observe its Behavior ---")

    # Simulation parameters
    num_steps = 2000
    start_position = 0
    
    # Store the path of the random walk
    path = [start_position]
    current_position = start_position
    
    for _ in range(num_steps):
        if random.random() < p_right:
            current_position += 1
        else:
            current_position -= 1
        path.append(current_position)

    print(f"Simulated {num_steps} steps of the random walk.")
    print("A plot of the walk's position over time will be generated.")
    print("For a non-positive recurrent chain, we expect it to wander far from the origin,")
    print("unlike a positive recurrent chain which would return to the origin frequently.")

    # Plotting the results
    plt.figure(figsize=(12, 6))
    plt.plot(path)
    plt.title("Symmetric Random Walk on Integers")
    plt.xlabel("Time Steps")
    plt.ylabel("Position")
    plt.grid(True)
    plt.show()

    print("\n--- Conclusion ---")
    print("The provided conditions are sufficient to conclude that the Markov chain is NOT positive recurrent.")


if __name__ == '__main__':
    demonstrate_non_positive_recurrence()
