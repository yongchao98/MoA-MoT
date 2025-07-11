import math

def count_cool_strings(n):
    """
    Calculates the number of 'cool strings' of maximal length 3n.
    The formula is derived from combinatorial arguments involving Standard Young Tableaux.
    Total = (2 * (3n)!) / ((n+1)! * (n+2)!)
    """
    if n <= 0:
        if n == 0:
            print("For n=0, there is 1 cool string (the empty string).")
            print("Formula: (2 * 0!) / (1! * 2!) = 1")
            return 1
        else:
            print("n must be a non-negative integer.")
            return 0

    try:
        # Calculate the components of the formula
        fact_3n = math.factorial(3 * n)
        fact_n_plus_1 = math.factorial(n + 1)
        fact_n_plus_2 = math.factorial(n + 2)

        # Calculate numerator and denominator
        numerator = 2 * fact_3n
        denominator = fact_n_plus_1 * fact_n_plus_2
        
        # The result must be an integer
        result = numerator // denominator

        # Print the detailed calculation as requested
        print(f"For n = {n}:")
        print(f"The number of strings is (2 * ({3*n})!) / (({n+1})! * ({n+2})!)")
        print(f"= (2 * {fact_3n}) / ({fact_n_plus_1} * {fact_n_plus_2})")
        print(f"= {numerator} / {denominator}")
        print(f"= {result}")
        return result

    except OverflowError:
        print(f"Calculation for n={n} resulted in numbers too large to compute.")
        return None
    except ValueError:
        print("Factorial is not defined for negative numbers.")
        return None


# Example usage for n=4
n_value = 4
count_cool_strings(n_value)

<<<11088>>>