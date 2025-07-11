import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for the described majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 0 or c <= 0:
        print("n and c must be positive integers.")
        return

    print(f"Calculating the mistake bound for n = {n} experts and c = {c} mistakes for removal.")

    if n < 3:
        # If n < 3, a wrong majority cannot be formed against the correct true expert.
        # Thus, mistakes of type M2 are not possible.
        # The bound is limited by the true expert's mistakes.
        bound = c - 1
        print("Since n < 3, a wrong majority cannot form against the correct true expert.")
        print(f"The bound is determined by the maximum mistakes of the true expert: M <= c - 1")
        print(f"Bound = {c} - 1")
        print(f"Upper Bound on Mistakes: {bound}")

    else:
        # For n >= 3, the bound is M <= M1 + M2
        # M1 <= c - 1
        # M2 <= floor((n - 1) * c / 2)
        m1_bound = c - 1
        m2_bound_val = (n - 1) * c / 2
        m2_bound = math.floor(m2_bound_val)
        total_bound = m1_bound + m2_bound

        print("The total bound M is the sum of two types of mistakes:")
        print("M1: Algorithm is wrong, and the true expert is wrong (M1 <= c - 1)")
        print("M2: Algorithm is wrong, but the true expert is right (M2 <= floor((n-1)*c / 2))")
        print("\nCalculating the bound:")
        print(f"Bound = (c - 1) + floor((n - 1) * c / 2)")
        print(f"Bound = ({c} - 1) + floor(({n} - 1) * {c} / 2)")
        print(f"Bound = {m1_bound} + floor({n-1} * {c} / 2)")
        print(f"Bound = {m1_bound} + floor({(n-1)*c} / 2)")
        print(f"Bound = {m1_bound} + floor({m2_bound_val})")
        print(f"Bound = {m1_bound} + {m2_bound}")
        print(f"Upper Bound on Mistakes: {total_bound}")


# Example Usage with placeholder values.
# You can change these values to see the bound for different scenarios.
n_experts = 10
c_mistakes = 5
calculate_mistake_bound(n_experts, c_mistakes)

final_n = n_experts
final_c = c_mistakes
if final_n < 3:
    final_bound = final_c - 1
else:
    final_bound = (final_c - 1) + math.floor(((final_n - 1) * final_c) / 2)
print(f'<<<{final_bound}>>>')