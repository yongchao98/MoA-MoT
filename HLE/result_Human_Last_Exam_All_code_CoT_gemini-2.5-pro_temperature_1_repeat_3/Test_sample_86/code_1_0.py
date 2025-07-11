import math

def calculate_min_prototypes():
    """
    Asks the user for N, then calculates the minimum total number of prototypes
    for a 1-NN classifier to perfectly classify N+1 concentric circles with
    balanced classes.
    """
    try:
        n_str = input("Enter the value for N (number of outer circles, must be an integer >= 0): ")
        N = int(n_str)
        if N < 0:
            print("Error: N must be a non-negative integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    print("\n--------------------")
    print("Step-by-step Calculation:")
    print("--------------------")

    if N == 0:
        print("For N=0, there is only one class (a single point at the origin).")
        print("To classify this single class, we only need one prototype.")
        total_prototypes = 1
        print("\nFinal Equation:")
        print(f"Total Prototypes = {total_prototypes}")
        
    else: # N >= 1
        # Step 1: Calculate k (prototypes per class)
        print("1. Determine the number of prototypes required per class, 'k'.")
        print("   To ensure perfect separation, k must be large enough to satisfy the condition")
        print("   derived from the geometry of the outermost circles.")
        print("   The minimum integer k is given by ceil(pi * N).")
        k = math.ceil(math.pi * N)
        print(f"   k = ceil(pi * {N}) = ceil({math.pi * N:.4f}) = {k}\n")

        # Step 2: Calculate number of classes
        print("2. Determine the total number of classes.")
        print("   The classes correspond to circles 0, 1, ..., N.")
        num_classes = N + 1
        print(f"   Number of classes = N + 1 = {N} + 1 = {num_classes}\n")

        # Step 3: Calculate total prototypes
        print("3. Calculate the total number of prototypes.")
        print("   This is the number of classes multiplied by the number of prototypes per class (k),")
        print("   as the problem requires balanced classes.")
        total_prototypes = num_classes * k
        print(f"   Total Prototypes = (Number of classes) * k = {num_classes} * {k} = {total_prototypes}\n")

        # Final result
        print("--------------------")
        print("Final Equation with values:")
        print("--------------------")
        print(f"Minimum Total Prototypes = ({N} + 1) * {k} = {total_prototypes}")

# Execute the function
calculate_min_prototypes()