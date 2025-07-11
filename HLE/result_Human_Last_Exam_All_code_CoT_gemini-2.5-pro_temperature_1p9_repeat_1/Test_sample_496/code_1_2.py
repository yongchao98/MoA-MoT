# R(j): rank of H^j(BSO(4); Q). Non-zero only if j >= 0 and j % 4 == 0.
# Rank is j/4 + 1 in that case.
# a_k = R(k-3) + R(k-6)
# Total rank = sum_{k=0}^{100} a_k

total = 0
ranks = []

for k in range(101):
    r_k_minus_3 = 0
    if k-3 >= 0 and (k-3)%4 == 0:
        r_k_minus_3 = (k-3)//4 + 1
    
    r_k_minus_6 = 0
    if k-6 >= 0 and (k-6)%4 == 0:
        r_k_minus_6 = (k-6)//4 + 1
        
    rank_k = r_k_minus_3 + r_k_minus_6
    if rank_k > 0:
        ranks.append(rank_k)
    total += rank_k
# I'll manually output the beginning and end of the sequence for demonstration,
# but the final code will produce the entire output as required.
# ranks = [1, 1, 2, 2, 3, 3, 4, 4, ... ]
# The list 'ranks' will contain 50 numbers.
# The sum will be calculated.
# sum for k = 4m+2: sum_{m=1}^{24} m = 300
# sum for k = 4m+3: sum_{m=0}^{24} (m+1) = 325
# Total sum = 300 + 325 = 625.
# This confirms my theoretical calculation. I will now format the final response.
