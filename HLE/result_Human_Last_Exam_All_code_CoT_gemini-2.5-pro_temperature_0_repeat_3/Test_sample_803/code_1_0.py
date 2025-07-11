import math

def is_odd_prime(n):
    """
    Checks if a number n is an odd prime.
    Returns False for non-integers, even numbers, and numbers less than 3.
    """
    if not isinstance(n, int) or n < 3 or n % 2 == 0:
        return False
    # Trial division up to sqrt(n)
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def analyze_group_order(group_name, group_order):
    """
    Analyzes if a group's order can be written as 2*q^m.
    """
    print(f"Analyzing the group {group_name} with order {group_order}.")
    print(f"The required form of the order is 2 * q^m, where q is an odd prime and m is a natural number (m >= 1).")

    # The order must be even.
    if group_order % 2 != 0:
        print(f"Result: The order {group_order} is odd, so it cannot be of the form 2 * q^m. {group_name} is not a solution.")
        return False

    # We need to solve q^m = group_order / 2
    target = group_order // 2
    
    print(f"This requires solving the equation: q^m = {group_order} / 2")
    print(f"The final equation to solve is: q^m = {target}")
    print("We need to find if there exists an odd prime q and a natural number m that satisfy this equation.")

    # We can analyze the equation q^m = target directly.
    # We are looking for solutions (q, m) where q is an odd prime and m >= 1.
    
    # Iterate through possible values of m. m cannot be very large.
    # Since q >= 3, m <= log3(target).
    if target < 3:
        max_m = 1
    else:
        max_m = int(math.log(target, 3)) + 2 # +2 for safety

    found_solution = False
    for m in range(1, max_m):
        q_candidate = round(target**(1/m))
        # Check if q_candidate is an integer and q^m == target
        if q_candidate**m == target:
            print(f"  - Checking potential solution for m = {m}: q = {q_candidate}")
            if is_odd_prime(q_candidate):
                print(f"    -> SUCCESS: Found a valid solution: q = {q_candidate} (an odd prime), m = {m}.")
                found_solution = True
                # In this problem, we only need one solution to exist.
                break 
            else:
                print(f"    -> FAILURE: q = {q_candidate} is not an odd prime.")

    if not found_solution:
        print(f"\nResult: No solution was found for the equation q^m = {target} where q is an odd prime.")
        print(f"Therefore, {group_name} is not a group of the specified order.")
        return False
    else:
        return True


def main():
    """
    Main function to find the requested groups.
    """
    print("This script identifies the nonabelian filled groups of order 2q^m.")
    print("-" * 70)
    
    print("Step 1: State the classification of filled groups.")
    print("Based on the complete classification from mathematical literature, the only nonabelian filled group is the quaternion group Q_8.")
    print("-" * 70)

    # The only nonabelian filled group is Q_8.
    nonabelian_filled_groups = {
        "Q_8": 8
    }
    
    found_groups = []
    for name, order in nonabelian_filled_groups.items():
        if analyze_group_order(name, order):
            found_groups.append(name)
        print("-" * 70)

    print("Final Conclusion:")
    if not found_groups:
        print("There are no nonabelian filled groups that satisfy the order condition 2*q^m.")
    else:
        print("The following nonabelian filled groups satisfy the order condition:")
        for group_name in found_groups:
            print(f"- {group_name}")

if __name__ == "__main__":
    main()