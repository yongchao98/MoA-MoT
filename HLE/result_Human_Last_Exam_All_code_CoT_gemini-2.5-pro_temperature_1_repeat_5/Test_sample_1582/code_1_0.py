import random
import matplotlib.pyplot as plt

def demonstrate_non_positive_recurrence():
    """
    Illustrates non-positive recurrence with a simple symmetric random walk.
    """
    # This Markov chain is a simple symmetric random walk on the integers (Z).
    # State space Sigma = Z. Transition p(i, i+1) = p(i, i-1) = 0.5.
    # This chain is irreducible and known to be null recurrent, so it is NOT
    # positive recurrent.

    # We can define a function f and a set A that satisfy the problem's conditions.
    # Let A = {-1, 0, 1}. A is a finite set.
    # Let f(x) = x^2. This function is non-negative and f(x) -> infinity as |x| -> infinity.

    # Let's check the submartingale condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A.
    # E[f(X_1) | X_0=x] = 0.5 * f(x+1) + 0.5 * f(x-1)
    #                    = 0.5 * (x+1)^2 + 0.5 * (x-1)^2
    #                    = 0.5 * (x^2 + 2x + 1) + 0.5 * (x^2 - 2x + 1)
    #                    = 0.5 * (2x^2 + 2) = x^2 + 1.
    # The drift is (x^2 + 1) - x^2 = 1, which is >= 0 for all x.
    # So the conditions of the problem are met.

    print("Demonstrating the submartingale property for a specific state:")
    x = 5
    f_x = x**2
    # The final equation is E[f(X_1)|X_0=x] = result
    expected_f_x_next = 0.5 * (x + 1)**2 + 0.5 * (x - 1)**2
    print(f"For state x = {x}, f(x) = {x**2} = {f_x}")
    print(f"E[f(X_1)|X_0={x}] = 0.5 * f({x+1}) + 0.5 * f({x-1}) = 0.5 * {(x+1)**2} + 0.5 * {(x-1)**2} = {expected_f_x_next}")
    drift = expected_f_x_next - f_x
    print(f"The drift is E[f(X_1)|X_0={x}] - f({x}) = {expected_f_x_next} - {f_x} = {drift}")
    print("The drift is non-negative, as required.\n")

    print("Simulating the random walk to observe its behavior.")
    print("A positive recurrent chain would return to the origin with a finite mean time.")
    print("A null recurrent or transient chain will wander off, and return times will tend to grow.")

    num_steps = 100000
    current_pos = 0
    path = [current_pos]
    return_times = []
    last_return_step = 0

    for i in range(1, num_steps + 1):
        if random.random() < 0.5:
            current_pos += 1
        else:
            current_pos -= 1
        path.append(current_pos)

        if current_pos == 0:
            return_time = i - last_return_step
            return_times.append(return_time)
            last_return_step = i
            
    print(f"\nSimulation finished after {num_steps} steps.")
    print(f"Number of returns to origin: {len(return_times)}")
    if return_times:
        print("Successive return times to the origin (excursion durations):")
        # To avoid flooding the output, we print the first 10 and last 10.
        if len(return_times) <= 20:
            print(return_times)
        else:
            print(f"First 10: {return_times[:10]}")
            print(f"Last 10: {return_times[-10:]}")
        print("\nNotice how the return times are not small and bounded, suggesting they might not have a finite average.")
    else:
        print("The chain did not return to the origin in the simulation time.")

if __name__ == '__main__':
    demonstrate_non_positive_recurrence()
