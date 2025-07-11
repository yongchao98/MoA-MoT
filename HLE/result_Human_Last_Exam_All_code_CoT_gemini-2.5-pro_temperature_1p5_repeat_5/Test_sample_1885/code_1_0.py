import random
from collections import Counter

def illustrate_bounding_principle():
    """
    This function illustrates a key principle used in proving that
    an uncountable subset of the functions must be bounded.
    """

    # Let's use finite numbers to represent the infinite cardinals
    # for the purpose of simulation. omega_2 is a larger cardinal than omega_1.
    OMEGA_2_SIZE = 1000  # Represents the number of functions, |omega_2|
    OMEGA_1_SIZE = 100   # Represents the number of possible values, |omega_1|

    print(f"Simulating with |ω₂| = {OMEGA_2_SIZE} functions and |ω₁| = {OMEGA_1_SIZE} possible values.\n")

    # We don't need to construct the full functions. We only need their values
    # at a specific coordinate, gamma. Let's pick one.
    gamma_coordinate = 42
    print(f"Let's examine the functions' values at a specific coordinate, γ = {gamma_coordinate}.\n")

    # For each function f_beta (where beta < omega_2), f_beta(gamma) is a value
    # less than omega_1. Let's generate these values randomly.
    # The actual sequence f_alpha has special properties, but for illustrating
    # this specific point, random values are sufficient.
    
    # values will store [f_0(gamma), f_1(gamma), ..., f_{OMEGA_2_SIZE-1}(gamma)]
    values_at_gamma = [random.randrange(OMEGA_1_SIZE) for _ in range(OMEGA_2_SIZE)]

    # According to the pigeonhole principle, since there are more functions than
    # possible values, some values must be repeated. Let's count them.
    value_counts = Counter(values_at_gamma)

    # Find the most common value and how many times it appears.
    most_common_value, times_appeared = value_counts.most_common(1)[0]
    
    print(f"Analyzing the values f_β({gamma_coordinate}) for all {OMEGA_2_SIZE} functions (β < {OMEGA_2_SIZE}):")
    print("-" * 50)
    print(f"The value that appears most often is: {most_common_value}")
    print(f"It appears {times_appeared} times.")
    
    # The theory for infinite cardinals guarantees that one value is taken
    # by |omega_2| many functions. Our simulation with finite numbers shows that
    # the number of repetitions is significantly larger than 1.
    # The expected number of repetitions for any value is OMEGA_2_SIZE / OMEGA_1_SIZE.
    expected_repetitions = OMEGA_2_SIZE / OMEGA_1_SIZE
    print(f"The guaranteed minimum number of repetitions for the most common value is ceil({OMEGA_2_SIZE}/{OMEGA_1_SIZE}) = {int(expected_repetitions + 0.999)}.")

    print("\nThis illustrates the core combinatorial argument:")
    print("For any coordinate γ, we can find a very large (uncountable, in the real problem) subset of functions")
    print("that all have the same value at γ. A full proof would weave these subsets together for all γ.")
    print("The existence of the bounding function g then follows.")


if __name__ == '__main__':
    illustrate_bounding_principle()
