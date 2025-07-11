import math

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a 
    federated knowledge graph system based on given parameters.
    """
    # Given parameters
    c = 0.95  # Identifier resolution consistency level
    b = 3     # Semantic version control branching factor

    # --- Model for FAIR components ---
    # The maximum theoretical score for each component is modeled as follows:

    # 1. Findability (f) is limited by the identifier consistency.
    f_max = c

    # 2. Accessibility (a) is assumed to be perfect in a best-case scenario.
    a_max = 1.0

    # 3. Interoperability (i) is limited by identifier consistency (c) and 
    #    semantic fragmentation (inversely related to branching factor b).
    i_max = c / b

    # 4. Reusability (r) is limited by the same factors as interoperability
    #    as it requires data to be findable and integrable.
    r_max = c / b

    # --- Calculate the final FAIR score R ---
    # R is the average of the four component scores.
    R_max = (f_max + a_max + i_max + r_max) / 4

    # --- Output the results ---
    print("Calculating the theoretical maximum FAIR compliance score (R):")
    print("-" * 60)
    print("The score (R) is modeled as the average of its four components:")
    print(f"R = (Findability + Accessibility + Interoperability + Reusability) / 4")

    print("\nComponent score models and values:")
    print(f"f = c = {f_max}")
    print(f"a = 1.0 (ideal implementation)")
    print(f"i = c / b = {c} / {b} = {i_max}")
    print(f"r = c / b = {c} / {b} = {r_max}")

    print("\nFinal equation with calculated values:")
    # The final equation as requested, showing all numbers.
    print(f"R = ({f_max} + {a_max} + {i_max} + {r_max}) / 4 = {R_max}")

if __name__ == "__main__":
    calculate_fair_score()
    
# The calculated value of R is (0.95 + 1.0 + 0.95/3 + 0.95/3) / 4 = 0.6458333333333333
# We will present this final value in the required format.
# <<<0.6458333333333333>>>