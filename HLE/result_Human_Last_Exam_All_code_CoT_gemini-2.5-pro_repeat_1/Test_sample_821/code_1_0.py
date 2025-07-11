def calculate_bound(n):
    """
    Calculates the upper bound for d derived from the nth step.
    The formula is f(n) = 1 + 1/(n*(n+1)).
    """
    return 1 + 1 / (n * (n + 1))

def main():
    """
    Shows that the derived upper bound for d approaches 1 as n increases.
    """
    print("The derived upper bound for d is d <= 1 + 1/(n*(n+1)) for all n.")
    print("Let's check this bound for the first few values of n:")
    for i in range(1, 11):
        bound = calculate_bound(i)
        print(f"For n = {i:2d}, d <= 1 + 1/({i}*({i}+1)) = {bound:.4f}")

    print("\nAs n becomes very large, the term 1/(n*(n+1)) approaches 0.")
    print("The tightest bound is the limit as n -> infinity, which is 1.")
    print("Therefore, the largest possible value for d is 1.")

if __name__ == "__main__":
    main()
