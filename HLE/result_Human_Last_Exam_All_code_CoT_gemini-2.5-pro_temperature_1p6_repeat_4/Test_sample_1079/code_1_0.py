import math

def main():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    
    # Properties of the H3 reflection group
    group_order = 120
    coxeter_number = 10
    degrees = [2, 6, 10]

    print("Step 1: The problem is to count elements `w` in H3 with `w*v = lambda*v`, where `lambda` is a primitive 10th root of unity and `v` is a regular eigenvector.")
    print(f"For H3, the Coxeter number is {coxeter_number}. A key theorem states that any eigenvector for an eigenvalue of order 10 is automatically regular.")
    print("Thus, the problem simplifies to counting elements that have a primitive 10th root of unity as an eigenvalue.\n")

    print("Step 2: Identify the relevant eigenvalues and apply Springer's theorem.")
    print(f"A primitive 10th root of unity is exp(2*pi*i*k/10) for k in {{1, 3, 7, 9}}.")
    print(f"The order of these eigenvalues is d={coxeter_number}. Since {coxeter_number} is a degree of H3 {degrees}, we can apply a theorem to count the elements.")
    print(f"The number of elements with a specific primitive {coxeter_number}th root of unity as an eigenvalue is |H3| / {coxeter_number}.\n")

    # The primitive roots can be grouped into complex conjugate pairs.
    # Pair 1: k=1, 9. Pair 2: k=3, 7.

    # Count for the first pair
    print("Step 3: Count the elements for the first pair of eigenvalues (k=1, 9).")
    num_elements_k1 = group_order // coxeter_number
    print(f"The number of elements with eigenvalue exp(2*pi*i*1/10) is {group_order} / {coxeter_number} = {num_elements_k1}.")
    print("These elements are distinct from those for other primitive roots.\n")
    
    # Count for the second pair
    print("Step 4: Count the elements for the second pair of eigenvalues (k=3, 7).")
    num_elements_k3 = group_order // coxeter_number
    print(f"The number of elements with eigenvalue exp(2*pi*i*3/10) is {group_order} / {coxeter_number} = {num_elements_k3}.")
    print("This set of elements is disjoint from the first set.\n")

    # Total count
    total_count = num_elements_k1 + num_elements_k3
    
    print("Step 5: The total number of such elements is the sum of the counts for the disjoint sets.")
    # The final print statement showing the equation as requested.
    print(f"Final Calculation: {num_elements_k1} + {num_elements_k3} = {total_count}")

if __name__ == "__main__":
    main()
