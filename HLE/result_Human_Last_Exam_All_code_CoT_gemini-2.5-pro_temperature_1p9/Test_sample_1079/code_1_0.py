import math

def calculate_elements():
    """
    Calculates the number of elements in the H3 reflection group with a regular
    eigenvector whose corresponding eigenvalue has order 10.
    """
    # Step 1: Define properties of the H3 reflection group.
    group_order = 120
    coxeter_number = 10

    # Step 2: Define a function to compute Euler's totient function phi(n).
    # This counts the number of conjugacy classes we are interested in.
    def euler_phi(n):
        """Computes Euler's totient function phi(n)."""
        count = 0
        for i in range(1, n):
            if math.gcd(i, n) == 1:
                count += 1
        return count

    # Step 3: Calculate the number of relevant conjugacy classes using phi(h).
    num_classes = euler_phi(coxeter_number)

    # Step 4: Calculate the size of each conjugacy class using |H3| / h.
    class_size = group_order // coxeter_number

    # Step 5: Calculate the total number of elements.
    total_elements = num_classes * class_size

    # Step 6: Print the explanation and results.
    print("The calculation follows the formula: phi(h) * (|H3| / h)")
    print(f"where |H3| = {group_order} is the group order and h = {coxeter_number} is the Coxeter number.")
    print("-" * 20)
    print(f"Number of relevant classes, phi({coxeter_number}) = {num_classes}")
    print(f"Size of each class, {group_order} / {coxeter_number} = {class_size}")
    print("-" * 20)
    print(f"Total number of elements = {num_classes} * {class_size} = {total_elements}")
    print("-" * 20)
    
    # As requested, outputting each number in the final equation:
    # phi(h) * (|H3| / h) = Total
    print("The numbers in the final equation are:")
    print(num_classes)
    print(group_order)
    print(coxeter_number)
    print(total_elements)

# Execute the function to get the answer.
calculate_elements()