import numpy as np

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit in the ring
    R = F_2[x,y]/(y^2 + x^4*y + x + 1).

    A unit u = a(x) + b(x)y has norm N(u) = a(x)^2 + a(x)b(x)x^4 + b(x)^2(x+1) = 1.
    The degree of u is max(deg(a), deg(b)+1).

    This function searches for the unit with the minimum degree by iterating through
    possible degrees and all polynomials of corresponding degrees.
    """
    P = np.polynomial.Polynomial
    one = P([1])
    x4 = P([0, 0, 0, 0, 1])
    x_plus_1 = P([1, 1])

    max_degree_to_check = 10
    for u_deg in range(1, max_degree_to_check):
        # The degree of u is max(deg(a), deg(b)+1)
        # Iterate over possible degrees for a(x) and b(x)
        # da is degree of a, db is degree of b. -1 means polynomial is 0.
        for da in range(-1, u_deg + 1):
            for db in range(-1, u_deg):
                if max(da if da!=-1 else -1, (db+1) if db!=-1 else 0) != u_deg:
                    continue

                # Iterate through all polynomials for given degrees
                # Number of polynomials of degree d is 2^d (leading coeff is 1)
                # But here we check all polynomials up to degree d.
                # So it's 2^(d+1) options
                
                # Number of candidates for a. Degree at most da.
                num_a = 2**(da + 1) if da >= 0 else 1
                # Number of candidates for b. Degree at most db.
                num_b = 2**(db + 1) if db >= 0 else 1

                for i in range(num_a):
                    if da == -1:
                        a_coeffs = [0]
                    else:
                        a_coeffs = [int(c) for c in bin(i)[2:].zfill(da+1)]
                    
                    # Ensure exact degree if da > -1
                    if da > -1 and len(a_coeffs) > 0 and a_coeffs[-1] == 0 and da > 0:
                        continue # Not exact degree
                    a = P(a_coeffs) % 2

                    for j in range(num_b):
                        if db == -1:
                            b_coeffs = [0]
                        else:
                            b_coeffs = [int(c) for c in bin(j)[2:].zfill(db+1)]
                        
                        # Ensure exact degree if db > -1
                        if db > -1 and len(b_coeffs)>0 and b_coeffs[-1] == 0 and db > 0:
                           continue # Not exact degree
                        b = P(b_coeffs) % 2
                        
                        # Skip trivial unit u=1 (a=1, b=0)
                        if (a==one and b == P([0])):
                           continue

                        # Calculate the norm
                        a_sq = (a * a) % 2
                        b_sq = (b * b) % 2
                        term1 = a_sq
                        term2 = (a * b * x4) % 2
                        term3 = (b_sq * x_plus_1) % 2
                        
                        norm = (term1 + term2 + term3) % 2
                        
                        # Check if the norm is 1
                        if norm == one:
                            a_str = ""
                            if da >=0 :
                                for d in range(da, -1, -1):
                                    if a.coef[d] == 1:
                                        if a_str: a_str += " + "
                                        if d == 0: a_str += "1"
                                        elif d == 1: a_str += "x"
                                        else: a_str += f"x^{d}"
                            if not a_str: a_str="0"
                            
                            b_str = ""
                            if db >= 0:
                               for d in range(db, -1, -1):
                                   if b.coef[d] == 1:
                                       if b_str: b_str += " + "
                                       if d == 0: b_str += "1"
                                       elif d == 1: b_str += "x"
                                       else: b_str += f"x^{d}"
                            if not b_str: b_str="0"

                            print(f"Found a unit u = a(x) + b(x)y with least degree {u_deg}.")
                            print(f"a(x) = {a_str}")
                            print(f"b(x) = {b_str}")
                            return u_deg

if __name__ == '__main__':
    degree = find_least_degree_unit()
    if degree:
        print(f"The least degree of a unit u != 1 is {degree}.")
