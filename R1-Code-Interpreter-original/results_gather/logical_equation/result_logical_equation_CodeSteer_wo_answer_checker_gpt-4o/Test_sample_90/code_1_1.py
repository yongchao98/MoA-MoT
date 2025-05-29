def solve():
    # Given numbers
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Helper function to check constraints
    def is_valid(A, B, C, D, E, F, G, H, I, J, K, L, M):
        return (
            E + K == 81 and
            K > J and
            E == 1.5 * G and
            B + J == 5 and
            M - F == -6 and
            K == 4.5 * M and
            M > B and
            F == 3.2 * C and
            C - E == -31 and
            L == 1.4 * C and
            A + M == 38 and
            F + I == 112
        )
    
    # Backtracking function
    def backtrack(assigned):
        if len(assigned) == 13:
            A, B, C, D, E, F, G, H, I, J, K, L, M = assigned
            if is_valid(A, B, C, D, E, F, G, H, I, J, K, L, M):
                print(f"<<<[{A}, {B}, {C}, {D}, {E}, {F}, {G}, {H}, {I}, {J}, {K}, {L}, {M}]>>>")
                return True
            return False
        
        for num in numbers:
            if num not in assigned:
                if backtrack(assigned + [num]):
                    return True
        return False
    
    # Start backtracking with an empty assignment
    backtrack([])

solve()