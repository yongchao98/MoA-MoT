import math

def demonstrate_infinite_transient_set():
    """
    This function demonstrates that one can construct an infinite transient set
    for the Doob's h-transform of a 2D simple random walk.

    The logic is that a set A = {y_k} is transient if the expected number of
    visits is finite. This expectation is proportional to the sum of 1/log(||y_k||).
    We choose ||y_k|| = exp(exp(k)), which makes the sum convergent.
    """

    print("To prove that not all transient sets are finite, we construct an infinite set A")
    print("and show that the expected number of visits to it is finite.")
    print("\nThe expected number of visits is proportional to the sum S, where:")
    print("S = Sum_{k=1 to inf} 1 / log(||y_k||)")
    print("We choose a set where ||y_k|| is approximately exp(exp(k)).")
    print("This simplifies the equation for each term of the sum:")
    print("Term_k = 1 / log(exp(exp(k))) = 1 / exp(k)")
    print("-" * 60)

    partial_sum = 0.0
    num_terms = 15

    print(f"We will now calculate the terms of the sum S = 1/e^1 + 1/e^2 + ...")
    print(f"and compute the partial sum for the first {num_terms} terms.")
    print("-" * 60)

    full_equation_str = []
    for k in range(1, num_terms + 1):
        term = 1 / math.exp(k)
        partial_sum += term
        
        # Output each number in the equation for the term
        print(f"Term for k={k}: 1 / e^{k} = {term:.8f}")
        
        if k <= 4:
            full_equation_str.append(f"{term:.4f}")

    # The problem asks to output the numbers in the final equation.
    # Here we show the start of the sum's equation.
    print("\nSo the final equation for the sum starts as:")
    print(f"S = {' + '.join(full_equation_str)} + ...")

    print("-" * 60)
    print(f"The partial sum of the first {num_terms} terms is: S_{num_terms} = {partial_sum:.8f}")
    
    exact_value = 1 / (math.e - 1)
    print(f"The exact value of the infinite sum is 1/(e-1) = {exact_value:.8f}")
    print("\nSince the sum converges, the expected number of visits is finite.")
    print("This means our constructed infinite set is transient.")
    print("Therefore, the statement that every transient set must be finite is false.")

demonstrate_infinite_transient_set()