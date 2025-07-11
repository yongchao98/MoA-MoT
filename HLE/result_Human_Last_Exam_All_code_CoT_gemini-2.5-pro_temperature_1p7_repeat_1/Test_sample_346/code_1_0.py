import math

def main():
    """
    Solves for log_2(n) where n is a Church integer represented by a
    complex SKI combinator expression.
    """
    # === Explanation of the Logic ===
    #
    # The SKI expression represents a Church numeral 'n'. Our goal is to find its value.
    # The core principles we'll use are:
    # 1.  Structure of the Expression: The top-level expression is S(x)(y)(z).
    #     This reduces to `x(z)(y(z))`. Given x=I, it becomes `I(z)(y(z))`, which is `z(y(z))`.
    #     So, the number n is `z` applied to the result of `y` applied to `z`.
    # 2.  Church Numerals Application: Applying Church numeral 'm' to 'n', written `m(n)`,
    #     is equivalent to exponentiation: n^m.
    # 3.  Key Combinators:
    #     - I: The identity combinator, representing the Church numeral 1.
    #     - S(I)(I): This combination reduces to the Church numeral 2.
    #     - S(S(K(S))(K)): This is the successor function, 'SUCC', where SUCC(n) = n + 1.

    # === Step-by-step Calculation ===

    print("Analyzing the SKI expression: S(I)(y)(z)")
    print("This reduces to n = z(y(z))\n")

    # Step 1: Calculate the value of z
    # z = S(S(K(S))(K))(S(I)(I))
    # This is SUCC applied to 2.
    church_1 = 1
    church_2 = 2
    z = church_2 + 1
    print(f"Step 1: Calculate z")
    print(f"z = SUCC(2) = 2 + 1 = {z}\n")

    # Step 2: Define the function y(arg)
    # y = S(I)(S(I)(K(C))) where C = S(S(K(S))(K))(I)
    # The application y(arg) reduces to arg(arg(C)).
    # C is SUCC applied to I (Church numeral 1).
    c = church_1 + 1
    print(f"Step 2: Define y(arg)")
    print(f"y(arg) reduces to arg(arg(C))")
    print(f"C = SUCC(1) = 1 + 1 = {c}")
    print(f"So, y(arg) = arg(arg(2))\n")

    # Step 3: Calculate y(z)
    # We substitute z=3 into y. So, y(3) = 3(3(2)).
    # Using the rule m(n) = n^m:
    # Inner part: 3(2) -> 2^3 = 8
    # Outer part: 3(8) -> 8^3 = 512
    inner_val_base = c
    inner_val_exp = z
    inner_val = inner_val_base ** inner_val_exp

    y_z_base = inner_val
    y_z_exp = z
    y_z = y_z_base ** y_z_exp

    print(f"Step 3: Calculate y(z) where z = {z}")
    print(f"y({z}) = {z}({z}({c}))")
    print(f"Inner part: {z}({c}) => {c}^{z} = {inner_val}")
    print(f"Outer part: {z}({inner_val}) => {inner_val}^{z} = {y_z}\n")

    # Step 4: Calculate n = z(y(z))
    # n = 3(512), which is 512^3.
    # We know 512 = 2^9, so n = (2^9)^3 = 2^(9*3) = 2^27.
    n_base = y_z
    n_exp = z
    log2_n_base = int(math.log2(n_base))

    print(f"Step 4: Calculate n")
    print(f"n = z(y(z)) = {z}({y_z})")
    print(f"This corresponds to {y_z}^{z}, which is {y_z}^{n_exp}.")
    print(f"Since {y_z} = 2^{log2_n_base}, n = (2^{log2_n_base})^{n_exp} = 2^({log2_n_base} * {n_exp})\n")

    # Step 5: Calculate the final answer, log_2(n)
    final_log_value = log2_n_base * n_exp
    print("Final step: Calculate log_2(n)")
    print(f"log_2(n) = log_2(2^{final_log_value})")
    
    # Final equation format as requested
    print("The final equation is:")
    print(f"log_2(({n_base})^{n_exp}) = log_2((2^{log2_n_base})^{n_exp}) = {log2_n_base} * {n_exp} = {final_log_value}")

if __name__ == "__main__":
    main()