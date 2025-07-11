import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.
    This is equivalent to finding b_4(Y), the number of senior conjugacy classes.
    """
    
    # Conjugacy classes of the icosahedral group A_5.
    # We describe them by element order, size, and the character value in the
    # standard 3D irreducible representation.phi = (1+sqrt(5))/2
    phi = (1 + math.sqrt(5)) / 2
    
    classes = [
        {'name': '1', 'order': 1, 'size': 1, 'char': 3, 'desc': 'Identity element'},
        {'name': '(12)(34)', 'order': 2, 'size': 15, 'char': -1, 'desc': 'Rotation by pi'},
        {'name': '(123)', 'order': 3, 'size': 20, 'char': 0, 'desc': 'Rotation by 2pi/3'},
        {'name': '(12345)', 'order': 5, 'size': 12, 'char': phi, 'desc': 'Rotation by 2pi/5'},
        {'name': '(13524)', 'order': 5, 'size': 12, 'char': 1 - phi, 'desc': 'Rotation by 4pi/5'},
    ]

    print("Calculating the age for each conjugacy class of the icosahedral group A_5:")
    print("-" * 70)

    num_senior_classes = 0

    for c in classes:
        print(f"Class represented by elements of type {c['name']} (order {c['order']}):")
        
        # The eigenvalues are roots of unity. Since det(g)=1, their product is 1.
        # The sum of eigenvalues is the character value.
        
        thetas = []
        eigenvalue_desc = ""

        if c['order'] == 1:
            # Identity element
            # Eigenvalues: 1, 1, 1
            thetas = [0.0, 0.0, 0.0]
            eigenvalue_desc = "1, 1, 1"
        elif c['order'] == 2:
            # Eigenvalues are +1 or -1. Sum is -1, product is 1. -> 1, -1, -1
            # -1 = exp(i*pi) = exp(2*pi*i * 1/2)
            thetas = [0.0, 0.5, 0.5]
            eigenvalue_desc = "1, -1, -1"
        elif c['order'] == 3:
            # Eigenvalues are 1, w, w^2 where w = exp(2*pi*i/3). Sum is 0, product is 1.
            thetas = [0.0, 1.0/3.0, 2.0/3.0]
            eigenvalue_desc = "1, exp(2*pi*i/3), exp(4*pi*i/3)"
        elif c['order'] == 5:
            # For A_5, the eigenvalues of a rotation are 1, exp(i*alpha), exp(-i*alpha).
            # The character value determines which rotation angle.
            if c['char'] > 0: # Corresponds to 1 + 2*cos(2*pi/5)
                # Angle is 2pi/5. Eigenvalues are 1, exp(2pi*i/5), exp(-2pi*i/5)=exp(8pi*i/5)
                thetas = [0.0, 1.0/5.0, 4.0/5.0]
                eigenvalue_desc = "1, exp(2*pi*i/5), exp(8*pi*i/5)"
            else: # Corresponds to 1 + 2*cos(4*pi/5)
                # Angle is 4pi/5. Eigenvalues are 1, exp(4pi*i/5), exp(-4pi*i/5)=exp(6pi*i/5)
                thetas = [0.0, 2.0/5.0, 3.0/5.0]
                eigenvalue_desc = "1, exp(4*pi*i/5), exp(6*pi*i/5)"
        
        age = sum(thetas)
        
        # We build the equation string for output.
        theta_strings = [f"{t:.2f}".rstrip('0').rstrip('.') for t in thetas]
        equation = " + ".join(theta_strings) + f" = {age:.0f}"

        print(f"  Eigenvalues: {eigenvalue_desc}")
        print(f"  Theta values (in [0,1)): {theta_strings[0]}, {theta_strings[1]}, {theta_strings[2]}")
        print(f"  Age calculation: {equation}")

        if age == 2:
            num_senior_classes += 1
            print("  This is a SENIOR class.")
        elif age == 1:
            print("  This is a JUNIOR class.")
        else: # age == 0 for identity
            print("  This is the IDENTITY class.")

        print("-" * 70)

    print("\nThe rank of H^2_c(Y, Q) is equal to b_4(Y), the number of senior (age=2) classes.")
    print(f"The number of senior classes found is: {num_senior_classes}")
    
    return num_senior_classes

if __name__ == '__main__':
    rank = calculate_rank()
    # The final answer format as requested.
    # The final print must contain the answer between <<<>>>
    print(f"\nFinal Answer: The rank of H^2_c(Y, Q) is {rank}.")
    print(f"<<<{rank}>>>")
