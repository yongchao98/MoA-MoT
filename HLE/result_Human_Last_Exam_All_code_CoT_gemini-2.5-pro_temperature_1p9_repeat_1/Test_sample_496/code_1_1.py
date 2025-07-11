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
    It prints each non-zero rank and then the final total rank.
    """
    total_rank = 0
    non_zero_ranks = []
    
    for k in range(101):
        rank_k = R(k - 3) + R(k - 6)
        if rank_k > 0:
            non_zero_ranks.append(rank_k)
        total_rank += rank_k

    for rank in non_zero_ranks:
        print(rank)
        
    print(total_rank)

if __name__ == "__main__":
    main()