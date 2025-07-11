import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.

    This rank is equal to the number of conjugacy classes of A_5 with age 1.
    """
    
    # The group G is the alternating group A_5. It has 5 conjugacy classes.
    # The representation is a 3D irreducible representation.
    # We list the classes and determine the age for each.
    # Note: phi is the golden ratio (1+sqrt(5))/2.
    
    # We will represent the eigenvalues conceptually. 
    # For an element of order n, eigenvalues are nth roots of unity, e.g., exp(2*pi*i*k/n).
    # The corresponding theta_k value is k/n.
    
    classes_info = [
        {
            "name": "Identity",
            "representative": "e",
            "order": 1,
            "eigenvalues": "1, 1, 1",
            "thetas": [0, 0, 0],
        },
        {
            "name": "3-cycles",
            "representative": "(1 2 3)",
            "order": 3,
            "eigenvalues": "1, w, w^2 (where w=exp(2*pi*i/3))",
            "thetas": [0, 1/3, 2/3],
        },
        {
            "name": "Products of two 2-cycles",
            "representative": "(1 2)(3 4)",
            "order": 2,
            "eigenvalues": "1, -1, -1",
            "thetas": [0, 1/2, 1/2],
        },
        {
            "name": "5-cycles (class 1)",
            "representative": "(1 2 3 4 5)",
            "order": 5,
            "eigenvalues": "1, z, z^4 (where z=exp(2*pi*i/5))",
            "thetas": [0, 1/5, 4/5],
        },
        {
            "name": "5-cycles (class 2)",
            "representative": "(1 2 3 5 4)",
            "order": 5,
            "eigenvalues": "1, z^2, z^3 (where z=exp(2*pi*i/5))",
            "thetas": [0, 2/5, 3/5],
        },
    ]

    print("Step-by-step calculation of the age for each conjugacy class of A_5:")
    print("-" * 70)

    age_one_classes_count = 0
    final_equation_terms = []

    for info in classes_info:
        name = info["name"]
        thetas = info["thetas"]
        age = sum(thetas)
        
        # Check if the sum is an integer for cleaner output
        if age == int(age):
            age = int(age)
        
        # The equation for age calculation
        thetas_str = " + ".join(map(str, info["thetas"]))
        
        print(f"Class: {name}")
        print(f"  Representative Element Order: {info['order']}")
        print(f"  Eigenvalues: {info['eigenvalues']}")
        print(f"  Corresponding Thetas (θ_j): {info['thetas']}")
        print(f"  Age = {thetas_str} = {age}\n")
        
        if age == 1:
            age_one_classes_count += 1
            final_equation_terms.append("1")

    print("-" * 70)
    print("The rank of H^2_c(Y, Q) is the number of conjugacy classes with age 1.")
    equation_str = " + ".join(final_equation_terms)
    print(f"Final Calculation: {equation_str} = {age_one_classes_count}")
    
    # Return the final answer in the required format
    return age_one_classes_count

if __name__ == '__main__':
    final_answer = solve_cohomology_rank()
    print(f"\n<<<__{final_answer}__>>>")
