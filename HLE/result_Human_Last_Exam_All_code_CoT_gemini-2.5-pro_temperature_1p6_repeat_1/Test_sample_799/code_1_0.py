import math

def demonstrate_sidon_dimension_limit():
    """
    This function demonstrates that the Hausdorff dimension of a constructed
    Sidon set can get arbitrarily close to 1/2.

    The dimension of the constructed set is given by the formula:
    dim = log(q) / log(2 * q^2)
    where q is the size of the underlying integer Sidon set.

    By taking the limit as q -> infinity, we can find the maximum dimension
    achievable with this method.
    limit(q->inf) [log(q) / log(2*q^2)]
    = limit(q->inf) [log(q) / (log(2) + 2*log(q))]
    Dividing the numerator and denominator by log(q):
    = limit(q->inf) [1 / (log(2)/log(q) + 2)]
    As q -> inf, log(q) -> inf, so log(2)/log(q) -> 0.
    The limit is 1 / (0 + 2) = 1/2.

    The following code calculates the dimension for a very large 'q' to
    numerically verify this limit.
    """
    
    # Choose a very large number for q to approximate infinity
    q = 10**200

    # Calculate the dimension using the formula
    dimension = math.log(q) / math.log(2 * q**2)
    
    print("Demonstration of the limit:")
    print(f"For a large construction parameter q = {q:.1e},")
    print(f"The calculated Hausdorff dimension is: {dimension}")
    print("\nThis demonstrates the dimension approaching the theoretical maximum.")

    # The exact maximum value is a known mathematical result.
    numerator = 1
    denominator = 2
    result = numerator / denominator

    print("\nThe theoretical maximum Hausdorff dimension is a fixed value.")
    print(f"The final answer is the fraction: {numerator} / {denominator}")
    print(f"Which equals: {result}")


demonstrate_sidon_dimension_limit()
