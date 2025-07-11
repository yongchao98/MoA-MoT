def solve():
    """
    This function determines the truth value of six statements about parenthesis strings.
    
    The reasoning is as follows:
    
    Let's analyze the statements by considering two extreme cases for the structure of the parenthesis string.
    Case A: A wide and shallow string, e.g., S_k = "(()()()...())" with k inner pairs.
    In S_k, we have one outer pair x_0 with L(x_0) = 2k+2, D(x_0) = 2, and k inner pairs x_i with L(x_i) = 2, D(x_i) = 1.
    
    Case B: A narrow and deep string, e.g., S'_k = "((...()))" with k nested pairs.
    In S'_k, we have k pairs x_j (j=1..k) with L(x_j) = 2j, D(x_j) = j.
    
    A crucial observation is the value of the function f(y) at y=1, which corresponds to the contribution of pairs like "()" with D=1.
    
    1. sum(log L(x)) = O(sum(log D(x))) -> False
       For S_k, LHS = log(2k+2) + k*log(2) which grows like O(k).
       RHS = log(2) + k*log(1) = log(2) (since log(1)=0), which is O(1).
       A function growing with k is not O(1).
    
    2. sum(log(log L(x))) = O(sum(log(log D(x)))) -> False
       Similar to (1). For S_k, assuming loglog(2) is a positive constant, LHS grows like O(k).
       The RHS sum is over D(x)>1 because log(log(1)) is undefined. So RHS is the constant log(log(2)).
       The statement is false.
    
    3. sum(log^5 L(x)) = O(sum(log^5 D(x))) -> False
       Same logic as (1). For S_k, LHS grows like O(k) and RHS is constant.
    
    4. sum(2^sqrt(log L(x))) = O(sum(2^sqrt(log D(x)))) -> True
       The key difference is that for D=1, the term is 2^sqrt(log 1) = 2^0 = 1, not 0.
       For S_k, LHS = 2^sqrt(log(2k+2)) + k * 2^sqrt(log 2). This is asymptotically O(k).
       RHS = 2^sqrt(log 2) + k * 2^sqrt(log 1) = 2^sqrt(log 2) + k. This is also O(k).
       The ratio LHS/RHS is bounded. For S'_k, RHS grows faster than LHS. So the statement holds.
    
    5. sum(L(x)^0.1) = O(sum(D(x)^0.11)) -> True
       For D=1, the term is 1^0.11 = 1.
       For S_k, LHS = (2k+2)^0.1 + k*2^0.1, which is O(k).
       RHS = 2^0.11 + k*1^0.11, which is O(k). The ratio is bounded.
       For S'_k, LHS is sum((2j)^0.1) ~ O(k^1.1), while RHS is sum(j^0.11) ~ O(k^1.11). Since 1.11 > 1.1, LHS is O(RHS). The statement holds.
    
    6. sum(L(x)^(1/4)) = O(sum(D(x)^(1/2))) -> True
       Same logic as (5). The exponent on D (0.5) is larger than on L (0.25).
       For S_k, both LHS and RHS are O(k).
       For S'_k, LHS ~ O(k^1.25) and RHS ~ O(k^1.5). LHS is O(RHS). The statement holds.
    
    The final result is a string of T/F values.
    """
    
    # The reasoning leads to the following conclusions:
    # 1. False
    # 2. False
    # 3. False
    # 4. True
    # 5. True
    # 6. True
    
    answer = "FFFTTT"
    print(answer)

solve()
print("<<<FFFTTT>>>")