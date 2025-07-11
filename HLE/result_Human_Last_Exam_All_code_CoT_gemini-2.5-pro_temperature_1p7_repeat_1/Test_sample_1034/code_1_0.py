import math

def solve_triangle_problem():
    """
    Solves the ambiguous triangle problem based on the given parameters.
    """
    # Step 1: Decode S from "SNYWR" using a phone keypad mapping.
    # P,Q,R,S -> 7
    # M,N,O   -> 6
    # W,X,Y,Z -> 9
    s_val = 76997
    print(f"The value of S is determined to be: {s_val}")

    # Step 2: Define the angle and side length.
    # BC = a = s_val
    # Angle B in degrees
    angle_b_deg = 30
    
    # The condition for two possible triangles is h < AC < BC,
    # where h is the altitude from C.
    # h = BC * sin(B) = s_val * sin(30) = s_val * 0.5
    
    # Step 3: Calculate the bounds for the length of side AC.
    lower_bound_float = s_val / 2.0
    upper_bound_float = s_val
    
    print(f"For two triangles to exist, the length of AC must be in the range ({lower_bound_float}, {upper_bound_float}).")

    # Step 4: Find the number of integers in this range.
    # The smallest integer AC can be is floor(lower_bound_float) + 1
    # The largest integer AC can be is ceil(upper_bound_float) - 1
    
    first_integer_ac = math.floor(lower_bound_float) + 1
    last_integer_ac = math.ceil(upper_bound_float) - 1

    print(f"The integer values for AC must be in the range [{first_integer_ac}, {last_integer_ac}].")

    # Step 5: Compute the count of these integer values.
    # The count is last_integer - first_integer + 1.
    count = last_integer_ac - first_integer_ac + 1
    
    print("\nThe number of integer values for AC is computed by the final equation:")
    print(f"Count = {last_integer_ac} - {first_integer_ac} + 1 = {count}")
    
    return count

# Execute the function and store the final answer.
final_answer = solve_triangle_problem()

print(f"\n<<< {final_answer} >>>")