def R(j):
    """
    Computes the rank of H^j(BSO(4); Q).
    This is non-zero only if j is a non-negative multiple of 4.
    """
    if j >= 0 and j % 4 == 0:
        return j // 4 + 1
    return 0

def main():
    """
    Calculates the total rank of A = H_{SO(4)}^*(SO(4) \ X) for degrees up to 100.
    The rank in degree k, denoted a_k, is given by the formula a_k = R_{k-3} + R_{k-6}.
    """
    total_rank = 0
    # The final equation string to be built
    final_equation = []
    
    for k in range(101):
        rank_k = R(k - 3) + R(k - 6)
        if rank_k > 0:
            final_equation.append(str(rank_k))
        total_rank += rank_k

    # Although not requested, showing the equation construction for clarity.
    # We will just print the numbers that are summed up.
    print("The total rank is the sum of the ranks in each degree from 0 to 100.")
    print("The ranks for degrees k where the rank is non-zero are:")
    
    # We are asked to output each number in the final equation.
    # We will print the non-zero ranks that are summed up.
    non_zero_ranks = " + ".join(final_equation)
    # The prompt asks us to print *each number* in the equation. Let's do that explicitly.
    print(f"The non-zero ranks are: {non_zero_ranks}")
    
    print("\nThe individual non-zero ranks being summed are:")
    for num in final_equation:
        print(num)

    print(f"\nThe total rank is: {total_rank}")

if __name__ == "__main__":
    main()
