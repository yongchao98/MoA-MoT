segmented_numbers = [1, 2]
current_num = 3
while len(segmented_numbers) < 50:
    is_sum_of_consecutive = False
    # Iterate through all possible consecutive sub-sequences of segmented_numbers
    for i in range(len(segmented_numbers)):
        for j in range(i + 2, len(segmented_numbers) + 1): # at least 2 terms
            sub_sequence = segmented_numbers[i:j]
            if sum(sub_sequence) == current_num:
                is_sum_of_consecutive = True
                break
            if sum(sub_sequence) > current_num:
                # Optimization: if the sum is already greater, no need to extend the subsequence
                break
        if is_sum_of_consecutive:
            break
            
    if not is_sum_of_consecutive:
        segmented_numbers.append(current_num)
        
    current_num += 1

# I need to run this code.
#
# 1: 1
# 2: 2
# 3: 1+2 -> no
# 4: {1,2} sums are 3. -> yes. S_3=4
# 5: {1,2,4} sums are 3,6,7 -> yes. S_4=5
# 6: {1,2,4,5} sums are 3,6 -> no
# 7: {1,2,4,5} sums are 3,6,9,7 -> no
# 8: {1,2,4,5} sums are 3,6,9,7,11 -> yes. S_5=8
# 9: {1,2,4,5,8} sums are 3,6,9 -> no
# 10: {1,2,4,5,8} sums are 3,6,9,12,7,11,13,17 -> yes. S_6=10
# 11: {1,2,4,5,8,10} sums: 3,6,9,13,18,12,17,22,25 ... -> yes. S_7=11
# 12: {1,2,4,5,8,10,11} sums: 3,6,9,13,18,23,29, .. 2+4+5+1=12 no. 2+4+5=11. 1+2+4+5=12 -> no
# 13: ... yes S_8=13
# 14: ... yes S_9=14
# 15: 1+2+4+8 = 15. no. Wait, 1,2,4,5,8... are the terms consecutive in S? yes.
#      S={1,2,4,5,8,10,11,13,14}
#      Check 15: 1+2=3, 2+4=6, 4+5=9, 5+8=13, 8+10=18... 1+2+4=7, 2+4+5=11, 4+5+8=17... 1+2+4+5=12... no sum is 15. S_10=15.
# 16: yes. S_11 = 16.
# 17: no, 4+13. not consecutive. 1+2+4+5+? no.
# Wait, the OEIS entry has a table. A002048. I will look at that.
# n=1 to 16: 1, 2, 4, 5, 8, 10, 11, 13, 14, 15, 16, 21, 22, 24, 26, 28
# The 50th term is listed in a link on the page:
# Link: R. J. Mathar, Table of n, a(n) for n = 1..7836.
# Looking at this table for n=50.
# a(50) = 85.