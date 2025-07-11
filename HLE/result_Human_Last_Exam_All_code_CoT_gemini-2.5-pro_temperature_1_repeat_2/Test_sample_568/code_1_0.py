import math

def solve_mixing_problem(epsilon, delta):
    """
    Calculates the number of repetitions needed for the liquids in two cups
    to become "the same" based on the mixing parameters epsilon and delta.

    The final derived formula for the number of steps t is:
    - If delta = 0.5, t = 1
    - If delta != 0.5, t = ceil(log(epsilon) / log(|1 - 2*delta|))
    """
    print(f"Analyzing the mixing problem with parameters:\nepsilon = {epsilon}\ndelta = {delta}\n")

    print("The condition for the liquids to be 'the same' is when the difference")
    print("in their compositions is less than or equal to epsilon. This leads to the inequality:")
    print("    |1 - 2*delta|^t <= epsilon")
    print(f"Plugging in the given values:")
    print(f"    |1 - 2*{delta}|^t <= {epsilon}\n")

    # Handle the special case where delta = 0.5
    if delta == 0.5:
        t = 1
        print("Since delta is 0.5, the base of the exponent is 0.")
        print("The liquids mix perfectly in a single step.")
        print("\nFinal Answer:")
        print(f"The number of repetitions required is: 1")
        return

    # General case for delta != 0.5
    base = abs(1 - 2 * delta)
    numerator = math.log(epsilon)
    denominator = math.log(base)
    t_float = numerator / denominator
    t_int = math.ceil(t_float)

    print("To solve for t, we use logarithms, which gives the equation:")
    print("    t >= log(epsilon) / log(|1 - 2*delta|)")
    print("\nCalculation steps:")
    # Using f-strings to format the output nicely.
    print(f"1. Calculate the numerator: log(epsilon) = log({epsilon}) = {numerator:.4f}")
    print(f"2. Calculate the base: |1 - 2*delta| = |1 - 2*{delta}| = {base:.4f}")
    print(f"3. Calculate the denominator: log(base) = log({base:.4f}) = {denominator:.4f}")
    print(f"4. Divide numerator by denominator: t >= {numerator:.4f} / {denominator:.4f} = {t_float:.4f}")
    print(f"5. Find the smallest integer t (the ceiling): ceil({t_float:.4f}) = {t_int}")

    print("\nFinal Answer:")
    print(f"The number of repetitions required is: {t_int}")

# --- Parameters ---
# You can change these values to see the result for different scenarios.
epsilon = 0.01  # The maximum allowed difference in composition (e.g., 0.01 for 1%)
delta = 0.1     # The fraction of liquid exchanged in each step

solve_mixing_problem(epsilon, delta)
