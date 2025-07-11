import math

def calculate_final_value():
    """
    This function calculates the final value based on the large-n approximation derived in the plan.
    """
    # The value of n from the problem statement
    n = 1000000000

    # a_i are the first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    # Based on the large-n approximation, l(a) ≈ -n/a.
    # The required sum is floor(sum_{i=1 to 10} l(a_i)).
    # This approximates to floor(-n * sum_{i=1 to 10} 1/a_i).

    # Calculate the sum of the reciprocals of the primes
    sum_of_reciprocals = 0
    for p in primes:
        sum_of_reciprocals += 1/p

    # Calculate the total sum using the approximation
    total_sum = -n * sum_of_reciprocals

    # The final result is the floor of the total sum
    final_answer = math.floor(total_sum)

    # Print the final equation showing the calculation as requested
    print("floor(", end="")
    print(f"{int(-n)} * (", end="")
    for i, p in enumerate(primes):
        print(f"1/{p}", end="")
        if i < len(primes) - 1:
            print(" + ", end="")
    print(f")) = {final_answer}")

calculate_final_value()