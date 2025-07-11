import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups to become "the same".

    Args:
        epsilon (float): The maximum allowed difference in concentration (between 0 and 1).
        delta (float): The fraction of liquid moved in each step (between 0 and 1).
    """

    print(f"Parameters: epsilon = {epsilon}, delta = {delta}\n")

    # The state of the system can be described by x_t, the fraction of red liquid in cup A at step t.
    # The recurrence relation is: x_{t+1} = x_t * (1 - 2*delta) + delta.
    # With x_0 = 1, the solution is: x_t = 0.5 + 0.5 * (1 - 2*delta)^t.
    #
    # The condition for the liquids to be "the same" is when the difference in their
    # red liquid concentration is at most epsilon.
    # Let r_A be the red fraction in A and r_B be the red fraction in B.
    # |r_A - r_B| <= epsilon
    # |x_t - (1 - x_t)| <= epsilon
    # |2*x_t - 1| <= epsilon
    #
    # Substituting the formula for x_t:
    # |2 * (0.5 + 0.5 * (1 - 2*delta)^t) - 1| <= epsilon
    # |1 + (1 - 2*delta)^t - 1| <= epsilon
    # |(1 - 2*delta)^t| <= epsilon
    # |1 - 2*delta|^t <= epsilon
    #
    # To solve for t, we take the logarithm of both sides:
    # t * log(|1 - 2*delta|) <= log(epsilon)
    #
    # Since |1 - 2*delta| < 1 (for delta in (0, 1)), its log is negative.
    # Dividing by a negative number flips the inequality:
    # t >= log(epsilon) / log(|1 - 2*delta|)
    #
    # The number of repetitions t must be the smallest integer satisfying this.
    # So, t = ceil(log(epsilon) / log(|1 - 2*delta|))

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        # If delta is 0.5, the liquids mix perfectly in one step.
        # The difference becomes 0, and for any epsilon > 0, the condition 0 <= epsilon is met.
        t = 1
        print("For delta = 0.5, the liquids mix perfectly in one step.")
        print(f"The number of repetitions needed is {t}.")
        return

    # For the general case, calculate the result using the formula.
    try:
        # These are the numbers that go into the final equation for t
        base = abs(1 - 2 * delta)
        log_epsilon = math.log(epsilon)
        log_base = math.log(base)

        print("The equation to solve for the number of repetitions t is:")
        print(f"t >= log(epsilon) / log(|1 - 2*delta|)")
        print("\nPlugging in the numbers:")
        print(f"t >= log({epsilon}) / log(|1 - 2*{delta}|)")
        print(f"t >= {log_epsilon} / {log_base}")
        
        # We need the smallest integer t, so we use math.ceil()
        t_float = log_epsilon / log_base
        print(f"t >= {t_float}")
        
        t = math.ceil(t_float)
        
        print(f"\nThe smallest integer number of repetitions is {t}.")

    except ValueError as e:
        print(f"Error calculating result: {e}")
        print("Please ensure epsilon and delta are in the range (0, 1).")
    except ZeroDivisionError:
        print("Error: Division by zero. This occurs if |1 - 2*delta| = 1 (i.e., delta=0 or delta=1).")
        print("For delta=0, no mixing occurs. For delta=1, the liquids just swap places.")
        print("In these cases, the mixture will never become 'the same' unless epsilon >= 1.")


# You can change these parameters to see the result for different values.
# As per the problem, epsilon and delta are in (0, 1).
epsilon_val = 0.05
delta_val = 0.1

solve_mixing_problem(epsilon_val, delta_val)