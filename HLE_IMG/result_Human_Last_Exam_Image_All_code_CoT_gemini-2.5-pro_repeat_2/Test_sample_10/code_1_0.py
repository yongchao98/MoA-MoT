import sys

def solve_and_predict_barrier():
    """
    This function predicts the inversion barrier for a PAH based on a quadratic model
    derived from the barriers of smaller homologues.
    """
    # Data points provided by the user
    # n: number of corannulene-like units
    # b: inversion barrier in kcal/mol
    n1, b1 = 1, 10
    n2, b2 = 2, 49
    
    # The target molecule has 3 units
    n_target = 3

    # We assume a quadratic relationship: B = a*n^2 + c
    # This gives a system of two linear equations:
    # 1) b1 = a * n1^2 + c
    # 2) b2 = a * n2^2 + c
    
    # We can solve this system for 'a' and 'c'.
    # Subtracting equation (1) from (2):
    # b2 - b1 = a * (n2^2 - n1^2)
    # So, a = (b2 - b1) / (n2^2 - n1^2)
    try:
        a = (b2 - b1) / (n2**2 - n1**2)
    except ZeroDivisionError:
        print("Error: n1 and n2 cannot be the same.", file=sys.stderr)
        return

    # Now substitute 'a' back into equation (1) to find 'c':
    # c = b1 - a * n1^2
    c = b1 - a * (n1**2)

    # Now we have the model, let's predict the barrier for the target molecule
    predicted_barrier = a * n_target**2 + c
    
    # Round to the nearest integer as requested
    final_result = int(round(predicted_barrier))

    # Print the explanation and the final calculation step-by-step
    print("A quadratic model (Barrier = a*n^2 + c) fits the provided data.")
    print(f"Based on the data for n=1 (Barrier={b1} kcal/mol) and n=2 (Barrier={b2} kcal/mol), the coefficients are a={int(a)} and c={int(c)}.")
    print("The resulting equation is: Barrier = 13 * n^2 - 3")
    print("\nFor the target molecule, triacenaphtho[...]triphenylene, n = 3.")
    print("The prediction calculation is:")
    
    # Output the final equation with all numbers, as requested.
    print(f"Barrier = {int(a)} * {n_target}^2 - {abs(int(c))} = {int(a)} * {n_target**2} - {abs(int(c))} = {int(a) * n_target**2} - {abs(int(c))} = {final_result}")
    
    print(f"\nThe predicted inversion barrier is {final_result} kcal/mol.")


# Execute the function
solve_and_predict_barrier()