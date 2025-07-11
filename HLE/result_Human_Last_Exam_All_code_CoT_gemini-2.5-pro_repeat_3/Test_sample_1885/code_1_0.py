import random
import math

def illustrate_pigeonhole_principle():
    """
    This function illustrates a key step in the proof.

    The proof relies on the fact that if you map a set of a larger
    cardinality (omega_2) to a set of a smaller cardinality (omega_1),
    then at least one element in the smaller set must be the image of a
    very large number of elements from the larger set.

    This code simulates this with large finite numbers.
    """
    # For illustration, let's use finite integers to represent the
    # cardinalities of the uncountable sets omega_2 and omega_1.
    omega_2_cardinality = 100000
    omega_1_cardinality = 1000

    print("--- Illustration of the Pigeonhole Principle Step ---")
    print(f"Simulating a mapping from a set of size {omega_2_cardinality} (like omega_2) to a set of size {omega_1_cardinality} (like omega_1).")

    # This dictionary will store the mapping from each 'alpha' to its value 'f_alpha(gamma_0)'
    mapping = {}
    for alpha in range(omega_2_cardinality):
        # We simulate f_alpha(gamma_0) by picking a random integer in the range of omega_1.
        value = random.randint(0, omega_1_cardinality - 1)
        mapping[alpha] = value

    # Now, we analyze the "fibers" of this mapping. A fiber is the set of all
    # 'alpha' values that map to the same 'delta_0' value.
    fibers = {}
    for alpha, value in mapping.items():
        if value not in fibers:
            fibers[value] = []
        fibers[value].append(alpha)

    # Find the value with the largest fiber.
    max_fiber_size = 0
    value_with_max_fiber = -1
    for value, fiber_list in fibers.items():
        if len(fiber_list) > max_fiber_size:
            max_fiber_size = len(fiber_list)
            value_with_max_fiber = value

    # The pigeonhole principle guarantees a minimum size for the largest fiber.
    guaranteed_min_size = math.ceil(omega_2_cardinality / omega_1_cardinality)

    print("\n--- Results ---")
    print(f"Found a value delta_0 = {value_with_max_fiber} that is the result of f_alpha(gamma_0) for {max_fiber_size} different alphas.")
    print(f"The generalized pigeonhole principle guarantees a fiber of at least size {guaranteed_min_size}.")

    print("\nThis demonstrates that an enormous number of functions must agree on the coordinate gamma_0.")
    print("In the actual proof with transfinite cardinals, mapping omega_2 to omega_1 guarantees a fiber of size omega_2, which is uncountable.")
    print("This uncountable set of agreeing functions is then used to derive a contradiction.")

if __name__ == '__main__':
    illustrate_pigeonhole_principle()
