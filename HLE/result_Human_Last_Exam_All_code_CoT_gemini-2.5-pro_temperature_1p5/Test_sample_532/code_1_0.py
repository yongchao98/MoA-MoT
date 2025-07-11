import math

def is_power_of_two(n):
    """Checks if a positive integer is a power of two."""
    if n <= 0:
        return False
    # A number n is a power of 2 if it's > 0 and n & (n-1) is 0.
    return (n & (n - 1)) == 0

def classify_finite_nilpotent_group(order, exponent):
    """
    Classifies a finite nilpotent group as filled or not based on a known theorem.

    According to the theorem, a finite nilpotent group G is filled if and only if
    it is NOT a 2-group with |G| > 2 and exp(G) > 2.

    Args:
        order (int): The order of the group |G|.
        exponent (int): The exponent of the group exp(G).

    Returns:
        bool: True if the group is filled, False otherwise.
    """
    is_2_group = is_power_of_two(order)

    # Check the condition for a group to be NOT filled
    if is_2_group and order > 2 and exponent > 2:
        return False
    else:
        return True

def main():
    """
    Demonstrates the classification with several example groups.
    """
    print("This script classifies finite nilpotent groups as 'filled' or 'not filled' based on their order and exponent.")
    print("-" * 20)

    # Example 1: Cyclic group C_3
    # order=3, exponent=3. Not a 2-group, so it is filled.
    order, exponent = 3, 3
    print(f"Group: C_3 (order={order}, exponent={exponent})")
    print(f"Result: {'Filled' if classify_finite_nilpotent_group(order, exponent) else 'Not Filled'}")
    print("-" * 20)

    # Example 2: Elementary abelian 2-group C_2 x C_2
    # order=4, exponent=2. It is a 2-group with order > 2, but its exponent is not > 2. So it is filled.
    order, exponent = 4, 2
    print(f"Group: C_2 x C_2 (order={order}, exponent={exponent})")
    print(f"Result: {'Filled' if classify_finite_nilpotent_group(order, exponent) else 'Not Filled'}")
    print("-" * 20)

    # Example 3: Cyclic group C_4
    # order=4, exponent=4. It is a 2-group with order > 2 and exponent > 2. So it is NOT filled.
    order, exponent = 4, 4
    print(f"Group: C_4 (order={order}, exponent={exponent})")
    print(f"Result: {'Filled' if classify_finite_nilpotent_group(order, exponent) else 'Not Filled'}")
    print("-" * 20)
    
    # Example 4: Quaternion group Q_8
    # order=8, exponent=4. It is a 2-group with order > 2 and exponent > 2. So it is NOT filled.
    order, exponent = 8, 4
    print(f"Group: Q_8 (order={order}, exponent={exponent})")
    print(f"Result: {'Filled' if classify_finite_nilpotent_group(order, exponent) else 'Not Filled'}")
    print("-" * 20)

    # Example 5: Cyclic group C_6 = C_2 x C_3
    # order=6, exponent=6. It is not a 2-group. So it is filled.
    order, exponent = 6, 6
    print(f"Group: C_6 (order={order}, exponent={exponent})")
    print(f"Result: {'Filled' if classify_finite_nilpotent_group(order, exponent) else 'Not Filled'}")
    print("-" * 20)


if __name__ == "__main__":
    main()
